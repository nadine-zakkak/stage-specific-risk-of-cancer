-- Purpose: Prepare table for cancer incidence analysis, using rectal bleeding & change in bowel habit cohorts
-- Author:  Nadine Zakkak 
-- Adopted from Becky White

-- --------------------------------------------------------------------------------
	-- Section A) Create Symptom event file
-- --------------------------------------------------------------------------------
/*
Using a join to the CPRD patient-level file, and practice-level file, 
rectal bleeding/change in bowel habit records occurring at any date in the ‘clinical’ file were selected, 
between 2007 and 2017
*/

-- identify patients that have any change in bowel habit (4) and/or rectal bleeding event (14)
drop temporary table if exists nadine.e2_initial_sx;
create temporary table nadine.e2_initial_sx as (
select c.d_cprdclin_key, c.e_patid, c.eventdate, c.medcode, l.symptom_number
from 18_299_Lyratzopoulos_e2.cprd_clinical c
inner join 18_299_Lyratzopoulos.lookup_core_symptom l on c.medcode = l.medcode
where l.symptom_number in (4, 14)
); 

select count(*), count(distinct e_patid) from nadine.e2_initial_sx; 
select count(*), count(distinct e_patid), symptom_number 
from nadine.e2_initial_sx 
group by symptom_number; 

create index `e_patid` on nadine.e2_initial_sx (`e_patid`);

-- only select patients from symptom cohort 
drop table if exists nadine.e2_cohort_initial_sx;
create table nadine.e2_cohort_initial_sx as (
select c.*
from nadine.e2_initial_sx c
inner join 18_299_Lyratzopoulos_e2.cprd_cohort_file l on c.e_patid = l.e_patid
);

select count(*), count(distinct e_patid) from nadine.e2_cohort_initial_sx; -- 519,670 rows & 344,770 unique e_patid
select count(*), count(distinct e_patid), symptom_number 
from nadine.e2_cohort_initial_sx 
group by symptom_number; 

create index `e_patid` on nadine.e2_cohort_initial_sx (`e_patid`);

-- created 'cprd_case_file' frmo patient, practice and linkage_eligiblity tables in script "linking_tables.sql"
drop table if exists nadine.e2_cohort_allsx;
create table nadine.e2_cohort_allsx
(
select d.d_cprdclin_key as origin_id, d.e_patid, 
d.eventdate, d.medcode, d.symptom_number as eventtype, 
l.dob, l.uts, l.crd, l.lcd, l.tod, l.deathdate, l.age30_date, l.age100_date, l.crd_oneyear, 
l.imd2015_10, l.gender
from nadine.e2_cohort_initial_sx d
LEFT JOIN nadine.e2_cprd_case_file l ON  d.e_patid = l.e_patid
WHERE 
 d.eventdate >= makedate(2007, 1)
AND d.eventdate <= makedate(2017, 365)
)
;
CREATE INDEX `eventdate` ON nadine.e2_cohort_allsx (`eventdate`);

select * from nadine.e2_cohort_allsx limit 20;

select count(*), count(distinct e_patid) from nadine.e2_cohort_allsx;
select count(*), count(distinct e_patid), eventtype 
from nadine.e2_cohort_allsx 
group by eventtype; 

-- --------------------------------------------------------------------------------
-- Section B) Flag potential index presentations
-- --------------------------------------------------------------------------------
/*
Rectal Bleeding/CIBH records were flagged as ‘potentially eligible’ if they also occurred:

-after the latest of
-- date patient was registered to the practice for at least one year (CRD + 1 year)
-- date practice's uts
-before the earliest of
-- practice's last collection date
-- patient's transfer out date
-- patient's death 
-between 30y and 99y old (inclusive) 

*/

-- Take off safe mode to enable column updates
set sql_safe_updates=0;

ALTER TABLE nadine.e2_cohort_allsx
ADD potent_elig INT(1)  DEFAULT NULL,
ADD sx INT(1) DEFAULT 1; 

-- Populate valid symptom flag = 1 if symptom event date met criteria
-- after follow up start, before follow up end, with (one day of follow up between start & end dates??)
UPDATE nadine.e2_cohort_allsx d
SET d.potent_elig = 1
WHERE 
 d.eventdate >= d.uts
AND d.eventdate >= d.crd_oneyear
AND  d.eventdate >= d.age30_date
AND d.eventdate <= d.lcd
AND (d.eventdate <= d.tod or d.tod is null)
AND (d.eventdate <= d.deathdate or d.deathdate is null)
AND d.eventdate < d.age100_date
; -- potential eligible events only occurred in patients from symptom cohort groups 

set sql_safe_updates=1;

select *
from  nadine.e2_cohort_allsx
limit 10;

select count(*), count(distinct e_patid) 
from nadine.e2_cohort_allsx
where potent_elig=1; 
select count(*), count(distinct e_patid), eventtype 
from nadine.e2_cohort_allsx
where potent_elig=1
group by eventtype; 

-- --------------------------------------------------------------------------------------------
-- C) Create patient cohort with 'potentially eligible' rectal bleeding/CIBH event from 2007-14
-- --------------------------------------------------------------------------------
/*
patients with at least one ‘potentially eligible’ rectal bleeding/CIBH record 
were selected as the ‘potentially eligible’ study cohort. 
*/

drop table if exists  nadine.e2_cohort_allsx_potelig_cohort;
create table nadine.e2_cohort_allsx_potelig_cohort
select distinct(e_patid)
from nadine.e2_cohort_allsx d
where potent_elig = 1;

CREATE INDEX `e_patid` ON nadine.e2_cohort_allsx_potelig_cohort (`e_patid`);
-- -----------------------------------------------------------------------
-- D) Select colon & rectal cancers for this cohort
-- --------------------------------------------------------------------------------
/*
all colon (11) and rectal (12) cancer diagnoses were selected from the Cancer Registry file, 
regardless of the date of occurrence. 
for patients with potential eligible rectal bleeding/CIBH diagnosis
*/

drop table if exists  nadine.e2_cohort_allsx_crc;
create table nadine.e2_cohort_allsx_crc
(
select d.e_patid, d.e_cr_patid, d.e_cr_id as origin_id, d.diagnosisdatebest as eventdate, d.site_icd10_o2, d.final_route, d.stage_best,
l2.icd10_3dig, l2.icd10_4dig, cancer_flag, l2.cancer_site_desc, l2.cancer_site_number as eventtype
from 18_299_Lyratzopoulos_e2.cancer_registration_tumour d
inner join nadine.e2_cohort_allsx_potelig_cohort l on d.e_patid = l.e_patid
inner join 18_299_Lyratzopoulos.lookup_core_cancersite l2 
on d.site_icd10_o2 = l2.icd10_4dig
where l2.cancer_site_number = 11 or l2.cancer_site_number = 12
)
;
CREATE INDEX `e_patid` ON nadine.e2_cohort_allsx_crc (`e_patid`);

set sql_safe_updates=0;
ALTER TABLE nadine.e2_cohort_allsx_crc
ADD sx INT(1) DEFAULT 0
;
set sql_safe_updates=1;

select * from nadine.e2_cohort_allsx_crc limit 20;
select count(*), count(distinct e_patid) from nadine.e2_cohort_allsx_crc; 
select count(*), count(distinct e_patid), eventtype from nadine.e2_cohort_allsx_crc group by eventtype;
select min(eventdate), max(eventdate) from nadine.e2_cohort_allsx_crc; 

-- ----------------------------------------------------------
-- E) Select all cancers except colorectal cancer for this cohort
-- --------------------------------------------------------------------------------
/*
all cancers except colroectal diagnoses were selected from the Cancer Registry file, 
regardless of the date of occurrence. 
for patients chosen to be in the final cohort 
*/

drop table if exists  nadine.e2_cohort_allsx_othercancer;
create table nadine.e2_cohort_allsx_othercancer
(
select d.e_patid, d.e_cr_patid, d.e_cr_id as origin_id, d.diagnosisdatebest as eventdate, d.site_icd10_o2, d.final_route, d.stage_best,
l2.icd10_3dig, l2.icd10_4dig, cancer_flag, l2.cancer_site_desc, l2.cancer_site_number as eventtype
from 18_299_Lyratzopoulos_e2.cancer_registration_tumour d
inner join nadine.e2_cohort_allsx_potelig_cohort l on d.e_patid = l.e_patid
inner join 18_299_Lyratzopoulos.lookup_core_cancersite l2 
on d.site_icd10_o2 = l2.icd10_4dig
where l2.cancer_site_number not in (11, 12) and l2.icd10_3dig like "%C%"
)
;
CREATE INDEX `e_patid` ON nadine.e2_cohort_allsx_othercancer (`e_patid`);

set sql_safe_updates=0;
ALTER TABLE nadine.e2_cohort_allsx_othercancer
ADD sx INT(1) DEFAULT 0
;
set sql_safe_updates=1;

select * from nadine.e2_cohort_allsx_othercancer limit 10;
select count(*), count(distinct e_patid) from nadine.e2_cohort_allsx_othercancer; -- 37290
select count(*), c.site_icd10_o2 from nadine.e2_cohort_allsx_othercancer c group by c.site_icd10_o2 order by c.site_icd10_o2;
select min(eventdate), max(eventdate) from nadine.e2_cohort_allsx_othercancer; -- 1995-01-03 to 2018-12-31

-- ----------------------------------------------
-- -- -- -- -- -- Continue in R -- -- -- -- -- --
-- ----------------------------------------------

-- -----------------------------------------------
-- Check if all good from R 
select * from nadine.e2_cohort_allsx_crc_final_dataset limit 50;
select count(*), count(distinct e_patid) from nadine.e2_cohort_allsx_crc_final_dataset; 
select count(*), crc from nadine.e2_cohort_allsx_crc_final_dataset group by crc;  
select count(e_patid) as c, e_patid from nadine.e2_cohort_allsx_crc_final_dataset group by e_patid having c > 1; 





