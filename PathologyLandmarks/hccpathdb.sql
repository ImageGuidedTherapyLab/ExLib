-- mysql  --local-infile < ./hccpathdb.sql
-- select * from Metadata.HCCPathdb;
-- set db
use Metadata;
-- NS database data
DROP TABLE IF EXISTS  Metadata.HCCPathdb;
CREATE TABLE Metadata.HCCPathdb(
id    bigint(20) NOT NULL AUTO_INCREMENT,
Rat                         VARCHAR(64)  NOT  NULL ,
TimePoint                   VARCHAR(64)       NULL ,
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



-- load dti data
LOAD DATA LOCAL INFILE './datalocation/CorrelativePathFiles.csv'
INTO TABLE Metadata.HCCPathdb
FIELDS TERMINATED BY ','  ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(Rat                         ,TimePoint                   ,T2WeightedReference         ,T2starMapOxygen             ,T2statMapMedicalAir         ,T2starMapAbsoluteChange     ,T2statMapPercentageChange   ,DCE                         ,T1PreContrast               ,T1PostContrast              ,PathologyHE                 ,PathologyPimo               )
SET metadata =  JSON_OBJECT( "Rat",Rat ,
                         "Time",cast(REPLACE(TimePoint,"weeks","") as Decimal)) ;

DROP PROCEDURE IF EXISTS HCCPathDBList ;
DELIMITER //
CREATE PROCEDURE HCCPathDBList 
( )
BEGIN
    SET SESSION group_concat_max_len = 10000000;
    select  concat("T2WeightedReference      =", group_concat( hcc.T2WeightedReference       separator ' ')) T2WeightedReference        from Metadata.HCCPathdb hcc;
    select  concat("T2starMapOxygen          =", group_concat( hcc.T2starMapOxygen           separator ' ')) T2starMapOxygen            from Metadata.HCCPathdb hcc;
    select  concat("T2statMapMedicalAir      =", group_concat( hcc.T2statMapMedicalAir       separator ' ')) T2statMapMedicalAir        from Metadata.HCCPathdb hcc;
    select  concat("T2starMapAbsoluteChange  =", group_concat( hcc.T2starMapAbsoluteChange   separator ' ')) T2starMapAbsoluteChange    from Metadata.HCCPathdb hcc;
    select  concat("T2statMapPercentageChange=", group_concat( hcc.T2statMapPercentageChange separator ' ')) T2statMapPercentageChange  from Metadata.HCCPathdb hcc;
    select  concat("DCE                      =", group_concat( hcc.DCE                       separator ' ')) DCE                        from Metadata.HCCPathdb hcc;
    select  concat("T1PreContrast            =", group_concat( hcc.T1PreContrast             separator ' ')) T1PreContrast              from Metadata.HCCPathdb hcc;
    select  concat("T1PostContrast           =", group_concat( hcc.T1PostContrast            separator ' ')) T1PostContrast             from Metadata.HCCPathdb hcc;
    select  concat("PathologyHE              =", group_concat( hcc.PathologyHE               separator ' ')) PathologyHE                from Metadata.HCCPathdb hcc;        
    select  concat("PathologyPimo            =", group_concat( hcc.PathologyPimo             separator ' ')) PathologyPimo              from Metadata.HCCPathdb hcc;
    select  concat("UpdateTransform          =", group_concat( hcc.UpdateTransform           separator ' ')) UpdateTransform            from Metadata.HCCPathdb hcc;
END //
DELIMITER ;
-- show create procedure Metadata.HCCPathDBList;
-- call Metadata.HCCPathDBList();
-- mysql  -sNre "call Metadata.HCCPathDBList();"
