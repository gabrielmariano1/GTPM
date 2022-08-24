-- EXECUTE IN MANAGEMENT STUDIO AND CLICK SAVE AS... TO SAVE THE RESULTS IN A TXT FILE
-- OPEN THE FILE USING NOTEPAD ++, REMOVE THE HEADERS AND CONVERT TO UTF-8 (IS UTF-9 BOM)

USE unitydb_logs

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

declare @gama int = 65 -- minimum time duration in seconds
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
SELECT DeviceLog_Today.Id AS EventID, DeviceLog_Today.UserId AS CHID, LocalTS AS EventDateTime, 
		CASE 
			WHEN EntityId IN(284,289,295,299) THEN 'B-T PAV JUNTA B BLOQ ENTRA GERAL IN' 
			WHEN EntityId IN (285,290,296,300) THEN 'B-T PAV JUNTA B BLOQ ENTRA GERAL OUT'
			WHEN EntityId IN (255,261,265) THEN 'C-T PAV JUNTA C BLQ FUNC IN'
			WHEN EntityId IN (256,262,266) THEN 'C-T PAV JUNTA C BLQ FUNC OUT'
		ELSE
			EntityName
		END AS CP
INTO #Temp
FROM [unitydb_logs].[dbo].[DeviceLog_Today]
JOIN [unitydb].[dbo].[Cards] ON Cards.Id = DeviceLog_Today.CredId
where EventId = 100
--and sourceid < 11
and (IsTemp <> 1 or istemp is null)
and LocalTS  BETWEEN '2020-03-01 00:00:00.000' AND '2020-03-30 23:59:59.999'  
order by eventdatetime asc

select COUNT(*) from #Temp

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
