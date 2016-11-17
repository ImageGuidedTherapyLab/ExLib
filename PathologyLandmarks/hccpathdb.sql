-- mysql  --local-infile < ./hccpathdb.sql
-- select * from HCCPath.metadata;
-- set db
use HCCPath;
-- NS database data
DROP TABLE IF EXISTS  HCCPath.metadata;
CREATE TABLE HCCPath.metadata(
id    bigint(20) NOT NULL AUTO_INCREMENT,
Rat                         VARCHAR(64)  NOT  NULL ,
TimePoint                   VARCHAR(64)       NULL ,
DataDir                     VARCHAR(256) NOT  NULL COMMENT 'Root data direcrtory',
T2WeightedReference         VARCHAR(256) NOT  NULL ,
T2WeightedReferenceLM       VARCHAR(256) GENERATED ALWAYS AS (REPLACE(T2WeightedReference,"RefImg.hdr","RefImgLM.hdr")) STORED,
UpdateTransform             VARCHAR(256) GENERATED ALWAYS AS (REPLACE(T2WeightedReference,"T2wReference/RefImg.hdr","updatetransform")) STORED,
T2starMapOxygen             VARCHAR(256)      NULL ,
T2statMapMedicalAir         VARCHAR(256)      NULL ,
T2starMapAbsoluteChange     VARCHAR(256)      NULL ,
T2statMapPercentageChange   VARCHAR(256)      NULL ,
DCE                         VARCHAR(256)      NULL ,
T1PreContrast               VARCHAR(256)      NULL ,
T1PostContrast              VARCHAR(256)      NULL ,
PathologyHE                 VARCHAR(256) NOT  NULL ,
PathologyHELM               VARCHAR(256) GENERATED ALWAYS AS (REPLACE(PathologyHE,"PathHE.svs","PathHELM.nii.gz")) STORED,
PathologyPimo               VARCHAR(256) NOT  NULL ,
PathologyPimoLM             VARCHAR(256) GENERATED ALWAYS AS (REPLACE(PathologyPimo,"PathPIMO.svs","PathPIMOLM.nii.gz")) STORED,
`metadata` JSON                                   NULL ,
PRIMARY KEY (id),
KEY `UID1` (`Rat`)
);

DROP PROCEDURE IF EXISTS HCCPath.ResetLabelStats ;
DELIMITER //
CREATE PROCEDURE HCCPath.ResetLabelStats()
BEGIN
  DROP TABLE IF EXISTS  HCCPath.lstat;
  CREATE TABLE HCCPath.lstat  (
   InstanceUID        VARCHAR(255)  NOT NULL COMMENT 'studyuid *OR* seriesUID', 
   SegmentationID     VARCHAR(80)   NOT NULL,  -- UID for segmentation file -- FIXME -- SOPUID NOT WORTH IT???  SegmentationSOPUID VARCHAR(255)   NOT NULL,  
   FeatureID          VARCHAR(80)   NOT NULL,  -- UID for image feature     -- FIXME -- SOPUID NOT WORTH IT???  FeatureSOPUID      VARCHAR(255)   NOT NULL,  
   LabelID            INT           NOT NULL,  -- label id for LabelSOPUID statistics of FeatureSOPUID      
   Mean               REAL              NULL,
   StdD               REAL              NULL,
   Max                REAL              NULL,
   Min                REAL              NULL,
   Count              INT               NULL,
   Volume             REAL              NULL,
   ExtentX            INT               NULL,
   ExtentY            INT               NULL,
   ExtentZ            INT               NULL,
   PRIMARY KEY (InstanceUID,SegmentationID,FeatureID,LabelID) );
END //
DELIMITER ;
show create procedure HCCPath.ResetLabelStats;
-- call HCCPath.ResetLabelStats();
show create table HCCPath.lstat;

DROP PROCEDURE IF EXISTS HCCPath.ResetOverlapStats ;
DELIMITER //
CREATE PROCEDURE HCCPath.ResetOverlapStats()
BEGIN
  DROP TABLE IF EXISTS  HCCPath.overlap;
  CREATE TABLE HCCPath.overlap(
   InstanceUID        VARCHAR(255)  NOT NULL COMMENT 'studyuid *OR* seriesUID',  
   FirstImage         VARCHAR(80)   NOT NULL,  -- UID for  FirstImage  
   SecondImage        VARCHAR(80)   NOT NULL,  -- UID for  SecondImage 
   LabelID            INT           NOT NULL,  -- label id for LabelSOPUID statistics of FeatureSOPUID      
   SegmentationID     VARCHAR(80)   NOT NULL,  -- UID for segmentation file  to join with lstat
   -- output of c3d firstimage.nii.gz secondimage.nii.gz -overlap LabelID
   -- Computing overlap #1 and #2
   -- OVL: 6, 11703, 7362, 4648, 0.487595, 0.322397  
   MatchingFirst      int           DEFAULT NULL,     --   Matching voxels in first image:  11703
   MatchingSecond     int           DEFAULT NULL,     --   Matching voxels in second image: 7362
   SizeOverlap        int           DEFAULT NULL,     --   Size of overlap region:          4648
   DiceSimilarity     real          DEFAULT NULL,     --   Dice similarity coefficient:     0.487595
   IntersectionRatio  real          DEFAULT NULL,     --   Intersection / ratio:            0.322397
   PRIMARY KEY (InstanceUID,FirstImage,SecondImage,LabelID) );
END //
DELIMITER ;
show create procedure HCCPath.ResetOverlapStats;
-- call HCCPath.ResetOverlapStats();
show create table HCCPath.overlap;


-- load dti data
LOAD DATA LOCAL INFILE './datalocation/CorrelativePathFiles.csv'
INTO TABLE HCCPath.metadata
FIELDS TERMINATED BY ','  ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(DataDir, Rat                         ,TimePoint                   ,T2WeightedReference         ,T2starMapOxygen             ,T2statMapMedicalAir         ,T2starMapAbsoluteChange     ,T2statMapPercentageChange   ,DCE                         ,T1PreContrast               ,T1PostContrast              ,PathologyHE                 ,PathologyPimo               )
SET metadata =  JSON_OBJECT( "Rat",Rat ,
                         "Time",cast(REPLACE(TimePoint,"weeks","") as Decimal)) ;

DROP PROCEDURE IF EXISTS HCCPath.HCCPathDBList ;
DELIMITER //
CREATE PROCEDURE HCCPath.HCCPathDBList 
( )
BEGIN
    SET SESSION group_concat_max_len = 10000000;
    select  concat("DATADIR                  =", group_concat( distinct hcc.DataDir                       )) DataDir                    from HCCPath.metadata hcc;
    select  concat("T2WeightedReference      =", group_concat( hcc.T2WeightedReference       separator ' ')) T2WeightedReference        from HCCPath.metadata hcc;
    select  concat("T2starMapOxygen          =", group_concat( hcc.T2starMapOxygen           separator ' ')) T2starMapOxygen            from HCCPath.metadata hcc;
    select  concat("T2statMapMedicalAir      =", group_concat( hcc.T2statMapMedicalAir       separator ' ')) T2statMapMedicalAir        from HCCPath.metadata hcc;
    select  concat("T2starMapAbsoluteChange  =", group_concat( hcc.T2starMapAbsoluteChange   separator ' ')) T2starMapAbsoluteChange    from HCCPath.metadata hcc;
    select  concat("T2statMapPercentageChange=", group_concat( hcc.T2statMapPercentageChange separator ' ')) T2statMapPercentageChange  from HCCPath.metadata hcc;
    select  concat("DCE                      =", group_concat( hcc.DCE                       separator ' ')) DCE                        from HCCPath.metadata hcc;
    select  concat("T1PreContrast            =", group_concat( hcc.T1PreContrast             separator ' ')) T1PreContrast              from HCCPath.metadata hcc;
    select  concat("T1PostContrast           =", group_concat( hcc.T1PostContrast            separator ' ')) T1PostContrast             from HCCPath.metadata hcc;
    select  concat("PathologyHE              =", group_concat( hcc.PathologyHE               separator ' ')) PathologyHE                from HCCPath.metadata hcc;        
    select  concat("PathologyPimo            =", group_concat( hcc.PathologyPimo             separator ' ')) PathologyPimo              from HCCPath.metadata hcc;
    select  concat("UpdateTransform          =", group_concat( hcc.UpdateTransform           separator ' ')) UpdateTransform            from HCCPath.metadata hcc;
END //
DELIMITER ;
-- show create procedure HCCPath.HCCPathDBList;
-- call HCCPath.HCCPathDBList();
-- mysql  -sNre "call HCCPath.HCCPathDBList();"
