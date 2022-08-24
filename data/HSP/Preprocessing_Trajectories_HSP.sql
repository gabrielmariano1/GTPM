USE unitydb_logs

if exists (select * from sysobjects where name='Trajectories' and xtype='U')
DROP TABLE Trajectories

create table Trajectories (
	CHID int,
	TrajectoryID int,
	SequenceID int,
	CP VARCHAR(50),
	EventDateTime datetime,
	Alpha int,
	Clusters varchar(800),
	Swarm varchar(400)
)


IF OBJECT_ID('tempdb.dbo.#Temp', 'U') IS NOT NULL DROP TABLE #Temp; 


;WITH CTE_1 AS
(
SELECT DeviceLog_Today.UserId AS CHID, DeviceLog_Today.Id AS EventID, LocalTS AS EventDateTime, A.CP, A.TimeStampID, A.GroupID, B.CHID AS ClusterCHID
FROM DeviceLog_Today
JOIN TimeStamps A ON DeviceLog_Today.Id = A.EventID
JOIN TimeStamps B ON A.TimeStampID = B.TimeStampID AND A.GroupID = B.GroupID
WHERE 1 = 1
--AND Events.CHID in (34, 37, 742,842,848,1156,17667,19770,28215,34634,50464,59659,62011,65499)
AND DeviceLog_Today.LocalTS  BETWEEN '2020-03-01 00:00:00.000' AND '2020-03-31 23:59:59.999' ---- 
), CTE_2 AS
(
SELECT CHID, EventDateTime, CP, TimeStampID, GroupID, Clusters = 
    STUFF((SELECT ',' + CAST(ClusterCHID AS VARCHAR(100))
           FROM CTE_1 b 
           WHERE b.CHID = a.CHID AND b.TimeStampID = a.TimeStampID AND b.GroupID = a.GroupID
          FOR XML PATH('')), 1, 1, '')
FROM CTE_1 a
--where TimeStampID IN (SELECT TimestampID FROM TimeStamps WHERE CHID  = 66371)
GROUP BY  CHID, EventDateTime, CP, TimeStampID, GroupID
)

SELECT CHID, EventDateTime, CP, TimeStampID, Clusters = 
    STUFF((SELECT ';' + CAST(Clusters AS VARCHAR(800))
           FROM CTE_2 b 
           WHERE b.CHID = a.CHID AND b.TimeStampID = a.TimeStampID
          FOR XML PATH('')), 1, 1, '')
INTO #Temp
FROM CTE_2 a
where CHID is not null
--where TimeStampID IN (SELECT TimestampID FROM TimeStamps WHERE CHID  = 66371)
GROUP BY  CHID, EventDateTime, CP, TimeStampID
ORDER BY EventDateTime DESC;

declare @listOfCHIDs table (CHID int);

DECLARE @CHID INT,
		@PreviousCHID INT = 0,
		@EventDateTime DATETIME,
		@PreviousEventDateTime DATETIME = 0,
		@FirstEventDateTime DATETIME = 0,
		@CP VARCHAR(50),
		@TimeStampID INT,
		@Clusters VARCHAR(800)


DECLARE @TrajectoryID INT = 1,
		@SequenceID INT = 0,
		@Alpha INT,
		@CurrentCluster VARCHAR(400),
		@Swarm VARCHAR(400)

DECLARE @CurrentScore INT,
		@HighestScore INT,
		@SwarmID INT

		
select COUNT(*) from #Temp

DECLARE timestamps_cursor CURSOR FOR 
SELECT	CHID,
		EventDateTime,
		CP,
		TimeStampID,
		Clusters 
FROM #Temp 
order by CHID ASC, EventDateTime asc

OPEN timestamps_cursor  


FETCH NEXT FROM timestamps_cursor   
INTO  @CHID, @EventDateTime, @CP, @TimeStampID, @Clusters

WHILE @@FETCH_STATUS = 0  
BEGIN 

	IF @CHID <>  @PreviousCHID 
	BEGIN
		SET @TrajectoryID = 1
		SET @SequenceID = 1
		SET @FirstEventDateTime = @EventDateTime
		SET @Alpha = 0
	END

	ELSE 
	BEGIN 
		IF DATEDIFF(HH, @PreviousEventDateTime, @EventDateTime) > 11
		BEGIN
			SET @TrajectoryID += 1
			SET @SequenceID = 1
			SET @FirstEventDateTime = @EventDateTime
			SET @Alpha = 0
		END
		ELSE 
		BEGIN
		SET @SequenceID += 1
		SET @Alpha = DATEDIFF(SECOND, @FirstEventDateTime, @EventDateTime)
		END
	END


	DECLARE Clusters_cursor CURSOR FOR 
	SELECT VALUE FROM SPLIT(';', @Clusters)

	OPEN Clusters_cursor  

	FETCH NEXT FROM Clusters_cursor   
	INTO  @CurrentCluster

	WHILE @@FETCH_STATUS = 0  
	BEGIN 

		DECLARE Swarms_cursor CURSOR FOR 
		SELECT DISTINCT SwarmID FROM Swarms WHERE CHID = @CHID

		OPEN Swarms_cursor  

		FETCH NEXT FROM Swarms_cursor   
		INTO  @SwarmID

		WHILE @@FETCH_STATUS = 0  
		BEGIN 

			WITH CTE_ AS
			(
				SELECT VALUE AS CHID FROM SPLIT(',', @CurrentCluster) INTERSECT SELECT CHID FROM Swarms WHERE SwarmID = @SwarmID
			)
			SELECT @CurrentScore = COUNT(1) FROM CTE_

			IF @CurrentScore > @HighestScore
			BEGIN
				SET @HighestScore = @CurrentScore
				SELECT @Swarm = 
				STUFF((SELECT ',' + CAST(CHID AS VARCHAR(100))
				FROM Swarms b WHERE SwarmID = @SwarmID
				FOR XML PATH('')), 1, 1, '')
				FROM Swarms a WHERE SwarmID = @SwarmID
			END

			FETCH NEXT FROM Swarms_cursor   
			INTO  @SwarmID
		END 

		CLOSE Swarms_cursor;  
		DEALLOCATE Swarms_cursor;

		FETCH NEXT FROM Clusters_cursor   
		INTO  @CurrentCluster
	END 

	IF @HighestScore < 2 SET @Swarm = NULL
	INSERT INTO Trajectories VALUES (@CHID, @TrajectoryID, @SequenceID, @CP, @EventDateTime, @Alpha, @Clusters, @Swarm)
	SET @HighestScore = 0

	CLOSE Clusters_cursor;  
	DEALLOCATE Clusters_cursor; 

	SET @PreviousCHID = @CHID
	SET @PreviousEventDateTime = @EventDateTime

	

	FETCH NEXT FROM timestamps_cursor   
	INTO  @CHID, @EventDateTime, @CP, @TimeStampID, @Clusters
END   
CLOSE timestamps_cursor;  
DEALLOCATE timestamps_cursor; 





--------------------------


select count(*) from Trajectories

--select * from Trajectories
--where CHID = 39668
--order by CHID asc, TrajectoryID asc, SequenceID asc

--select distinct chid from Trajectories


select top 1000 * from #Temp
order by TimeStampID