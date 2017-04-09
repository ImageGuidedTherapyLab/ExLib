-- mysql  --local-infile < ./hccpathdb.sql
-- select * from HCCPath.metadata;
-- set db
use HCCPath;
-- NS database data
DROP TABLE IF EXISTS  HCCPath.metadata;
CREATE TABLE HCCPath.metadata(
id    bigint(20) NOT NULL AUTO_INCREMENT,
Rat                         VARCHAR(64)  NOT  NULL ,
Status                              INT       NULL COMMENT '0 = control, 1= treatment',
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
(DataDir, Rat, Status,TimePoint                   ,T2WeightedReference         ,T2starMapOxygen             ,T2statMapMedicalAir         ,T2starMapAbsoluteChange     ,T2statMapPercentageChange   ,DCE                         ,T1PreContrast               ,T1PostContrast              ,PathologyHE                 ,PathologyPimo               )
SET metadata =  JSON_OBJECT( "Rat",Rat ,
                         "Time",cast(REPLACE(TimePoint,"weeks","") as Decimal)) ;

DROP PROCEDURE IF EXISTS HCCPath.HCCPathDBList ;
DELIMITER //
CREATE PROCEDURE HCCPath.HCCPathDBList 
(IN csvfile  varchar(255))
BEGIN
-- FIXME replace csvfile  into  table
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


DROP PROCEDURE IF EXISTS HCCPath.HCCPathOutput ;
DELIMITER //
CREATE PROCEDURE HCCPath.HCCPathOutput 
(IN radius int)
BEGIN
-- set  @varRadius=20;   collate utf8_unicode_ci;
set  @varRadius=radius ;
select md.Rat, md.Status, md.TimePoint,
       md.PathologyHE, hee.mean as EntropyHE,heh.mean as HaralickHE, hev.mean as DistViableHE, hen.mean as DistNecrosisHE, 
       md.PathologyPimo ,pie.mean as EntropyPimo,pih.mean as HaralickPimo, pio.mean as DistO2Pimo,
       md.DCE ,dce.mean as DCEMean,avg.mean as DCEAvg,
       md.T2WeightedReference , t2abs.mean T2Abs, t2air.mean MedAir, t2oxy.mean OxyT2, t2pct.mean T2Pct,
       t1pre.mean T1Pre, t1pst.mean T1Pst
from        HCCPath.metadata md
left  join  HCCPath.lstat    heh    on (md.Rat=heh.InstanceUID   and heh.LabelID   = 4 and heh.FeatureID=CONCAT('PathHE000.HaralickCorrelation_',@varRadius ,'.nii.gz'))
left  join  HCCPath.lstat    pih    on (md.Rat=pih.InstanceUID   and pih.LabelID   = 4 and pih.FeatureID=CONCAT('PathPIMO000.HaralickCorrelation_',@varRadius ,'.nii.gz')) 
left  join  HCCPath.lstat    hev    on (md.Rat=hev.InstanceUID   and hev.LabelID   = 1 and hev.FeatureID=CONCAT('PathHELMdist.nii.gz'))
left  join  HCCPath.lstat    hen    on (md.Rat=hen.InstanceUID   and hen.LabelID   = 3 and hen.FeatureID=CONCAT('PathHELMdist.nii.gz'))
left  join  HCCPath.lstat    pio    on (md.Rat=pio.InstanceUID   and pio.LabelID   = 1 and pio.FeatureID=CONCAT('PathPIMOLMdist.nii.gz'))
left  join  HCCPath.lstat    hee    on (md.Rat=hee.InstanceUID   and hee.LabelID   = 4 and hee.FeatureID=CONCAT('PathHE000.Entropy_',@varRadius ,'.nii.gz'))
left  join  HCCPath.lstat    pie    on (md.Rat=pie.InstanceUID   and pie.LabelID   = 4 and pie.FeatureID=CONCAT('PathPIMO000.Entropy_',@varRadius ,'.nii.gz')) 
left  join  HCCPath.lstat    t2abs  on (md.Rat=t2abs.InstanceUID and t2abs.LabelID = 4 and t2abs.FeatureID=CONCAT('T2AbsChange.hdr')) 
left  join  HCCPath.lstat    t2air  on (md.Rat=t2air.InstanceUID and t2air.LabelID = 4 and t2air.FeatureID=CONCAT('MedAirT2.hdr')) 
left  join  HCCPath.lstat    t2oxy  on (md.Rat=t2oxy.InstanceUID and t2oxy.LabelID = 4 and t2oxy.FeatureID=CONCAT('OxyT2.hdr')) 
left  join  HCCPath.lstat    t2pct  on (md.Rat=t2pct.InstanceUID and t2pct.LabelID = 4 and t2pct.FeatureID=CONCAT('T2PctChange.hdr')) 
left  join  HCCPath.lstat    t1pre  on (md.Rat=t1pre.InstanceUID and t1pre.LabelID = 4 and t1pre.FeatureID=CONCAT('T1pre.hdr')) 
left  join  HCCPath.lstat    t1pst  on (md.Rat=t1pst.InstanceUID and t1pst.LabelID = 4 and t1pst.FeatureID=CONCAT('T1post.hdr')) 
left  join  HCCPath.lstat    avg    on (md.Rat=avg.InstanceUID   and avg.LabelID   = 4 and avg.FeatureID=CONCAT('DCEavg.hdr')) 
left  join  HCCPath.lstat    dce    on (md.Rat=dce.InstanceUID   and dce.LabelID   = 4 and dce.FeatureID=CONCAT('DCE.hdr'));
END //
DELIMITER ;
-- show create procedure HCCPath.HCCPathOutput;
-- call HCCPath.HCCPathOutput(20);
-- mysql  -re "call HCCPath.HCCPathOutput(20);" | sed "s/\t/,/g;s/NULL//g" > analysissummary.csv
