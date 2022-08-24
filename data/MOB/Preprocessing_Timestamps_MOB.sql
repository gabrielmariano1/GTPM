
if exists (select * from sysobjects where name='TimeStamps' and xtype='U')
DROP TABLE TimeStamps

create table TimeStamps (
	EventDateTime DateTime,
	CP VARCHAR(50),
	TimeStampID int, 
	GroupID int,
	CHID int,
	EventID bigint
)
go

CREATE NONCLUSTERED INDEX [IX_CHID_TIMESTAMP_GROUPID] ON [dbo].[TimeStamps]
(
	[CHID] ASC
)
INCLUDE ( 	[TimeStampID],
	[GroupID]) WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, SORT_IN_TEMPDB = OFF, DROP_EXISTING = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
GO

declare @gama int = 30 -- minimum time duration in seconds
declare @listOfCHIDs table (EventID int, CHID int);
declare @LastlistOfCHIDs table (EventID int, CHID int);
declare @timestamp int;

DECLARE @EventID BIGINT,
		@CHID INT,
		@EventDateTime DATETIME,
		@LastEventDateTime DATETIME,
		@CP VARCHAR(50),
		@LastCP VARCHAR(50),
		@TimeStampID INT = 0,
		@ClusterID INT = 0,
		@ID INT


IF OBJECT_ID('tempdb.dbo.#Temp', 'U') IS NOT NULL DROP TABLE #Temp; 
SELECT	EventID, CHID, DATEADD(MI, GMTOffSet, EventDateTime) AS EventDateTime, 
		(SELECT CASE SourceID
		WHEN 1 THEN 'Reception'
		WHEN 2 THEN 'Reception'
		WHEN 3 THEN 'Reception'
		WHEN 15 THEN 'Warehouse'
		WHEN 17 THEN 'SS'
		WHEN 19 THEN 'G1'
		WHEN 21 THEN 'G1'
		WHEN 23 THEN 'G2'
		WHEN 25 THEN 'G2'
		WHEN 27 THEN 'G3'
		WHEN 28 THEN 'G3'
		WHEN 31 THEN 'Administration'
		WHEN 33 THEN 'Emergency'
		WHEN 37 THEN 'Turnstile Entry'
		WHEN 38 THEN 'Turnstile Exit'
		WHEN 39 THEN 'Gate 1SS'
		END) + ' '  + ISNULL((SELECT CASE EventHWID
		WHEN 73 THEN 'Entry' --Entry
		WHEN 75  THEN 'Exit' --Exit
		END), '' ) AS CP
INTO #Temp
FROM [W_Access_Events_Passarelli].[dbo].[Events]
where sourcetype = 1 
--and sourceid < 11
and chtype = 2
and eventhwid IN (71, 73, 75)
and EventDateTime  BETWEEN '2018-03-01 00:00:00.000' AND '2018-03-30 23:59:59.999'   --IN ('2018/03/26', '2018/03/27', '2018/03/28', '2018/03/29', '2018/03/30')
--and cast(eventdatetime as date) BETWEEN '2018/03/01' AND '2018/03/30'   --IN ('2018/03/26', '2018/03/27', '2018/03/28', '2018/03/29', '2018/03/30')
--and eventhwid = 75
--and chid in (select chid from cte_chid)
--and eventdatetime > '2018-03-26 18:50:00.000'
order by eventdatetime asc

select * from #Temp

DECLARE events_cursor CURSOR FOR 
SELECT EventID, CHID, EventDateTime, CP FROM #Temp order by CP, eventdatetime asc
OPEN events_cursor  


FETCH NEXT FROM events_cursor   
INTO @EventID, @CHID, @EventDateTime, @CP

WHILE @@FETCH_STATUS = 0  
BEGIN

	--Set the list of chids of the group
	DELETE FROM @listOfCHIDs
	INSERT INTO @listOfCHIDs
	SELECT EventID, CHID
	FROM  #Temp
	WHERE CP = @CP
	AND EventDateTime BETWEEN  @EventDateTime AND DATEADD(SECOND, @gama, @EventDateTime)

	--SELECT CHID FROM @listOfCHIDs
	--SELECT CHID FROM @LastlistOfCHIDs
	--SELECT CHID FROM @listOfCHIDs EXCEPT SELECT CHID FROM @LastlistOfCHIDs
	--SELECT * FROM TimeStamps
	--SELECT '----'


	-- is this a new Group from the same Timestamp?
	IF @CP = @LastCP AND DATEDIFF(SECOND, @LastEventDateTime, @EventDateTime) < @gama AND EXISTS (SELECT EventID, CHID FROM @listOfCHIDs EXCEPT SELECT EventID, CHID FROM @LastlistOfCHIDs)
	BEGIN
		SET @ClusterID += 1
		INSERT INTO TimeStamps
		SELECT @EventDateTime, @CP, @TimeStampID, @ClusterID, CHID, EventID FROM @listOfCHIDs
	END
	
	IF (@CP <> @LastCP OR DATEDIFF(SECOND, @LastEventDateTime, @EventDateTime) > @gama) 
	BEGIN
		SET @ClusterID = 1
		SET @TimeStampID += 1
		INSERT INTO TimeStamps
		SELECT @EventDateTime, @CP, @TimeStampID, @ClusterID, CHID, EventID FROM @listOfCHIDs
	END

	


	SET @LastCP = @CP
	SET @LastEventDateTime = @EventDateTime

	DELETE FROM @LastlistOfCHIDs
	INSERT INTO @LastlistOfCHIDs SELECT EventID, CHID FROM  @listOfCHIDs
	

	FETCH NEXT FROM events_cursor   
	INTO @EventID, @CHID, @EventDateTime, @CP
END   
CLOSE events_cursor;  
DEALLOCATE events_cursor; 

--select * from #temp
----WHERE CHID = 24144
--order by eventdatetime asc

SELECT * FROM TimeStamps
order by TimeStampID asc


