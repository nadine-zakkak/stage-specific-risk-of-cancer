-- Purpose: Prepare table for cancer stage analysis, using colon and rectal cancer cohorts
-- Author:  Nadine Zakkak 

-- ------------------------------------------------------------------------------------------
-- Step 1: Find all colon and rectal cancers for all patients between 1/1/2012 and 31/12/2018
-- ------------------------------------------------------------------------------------------
drop table if exists nadine.e2_allcrc;
create table nadine.e2_allcrc
select  d.e_patid, d.e_cr_patid, d.e_cr_id as origin_id, d.diagnosisdatebest as eventdate, d.final_route, d.stage_best,
l2.icd10_3dig, l2.icd10_4dig, cancer_flag, l2.cancer_site_desc, l2.cancer_site_number as eventtype, 
l1.dob, l1.uts, l1.crd, l1.lcd, l1.tod, l1.deathdate, l1.age30_date, l1.age100_date, l1.crd_oneyear, 
l1.imd2015_10, l1.gender, l1.random_sample, l1.cohort -- !!!! cohort is from ANY symptom cohort)
from 18_299_Lyratzopoulos_e2.cancer_registration_tumour d
inner join nadine.e2_cprd_case_file l1 on d.e_patid = l1.e_patid
inner join 18_299_Lyratzopoulos.lookup_core_cancersite l2 on d.site_icd10_o2 = l2.icd10_4dig
where (l2.cancer_site_number = 11 or l2.cancer_site_number = 12)
and d.diagnosisdatebest >= makedate(2012, 1)
and d.diagnosisdatebest <= makedate(2018, 365)
;

select * from nadine.e2_allcrc limit 10;
select count(*), count(distinct e_patid) from nadine.e2_allcrc; 
select count(*), count(distinct e_patid), cancer_site_desc from nadine.e2_allcrc group by cancer_site_desc; 
-- only random
select count(*), count(distinct e_patid), cancer_site_desc from nadine.e2_allcrc where random_sample = 1 and cohort= 0 group by cancer_site_desc; 
-- only cohort
select count(*), count(distinct e_patid), cancer_site_desc from nadine.e2_allcrc where random_sample = 0 and cohort = 1 group by cancer_site_desc; 
-- both
select count(*), count(distinct e_patid), cancer_site_desc from nadine.e2_allcrc where random_sample = 1 and cohort = 1 group by cancer_site_desc; 
select min(eventdate), max(eventdate) from nadine.e2_allcrc; 

CREATE INDEX `e_patid` ON nadine.e2_allcrc (`e_patid`);

-- ------------------------------------------------------------------------------------------
-- Step 2: Find all rectal bleeding and change in bowel habit symptom for identified patients with CRC
-- ------------------------------------------------------------------------------------------
-- identify patients that have any change in bowel habit (4) and/or rectal bleeding event (14)
drop table if exists nadine.e2_allcrc_initial_sx_new;
create table nadine.e2_allcrc_initial_sx_new as (
select distinct c.d_cprdclin_key as origin_id, c.e_patid, c.eventdate, l.symptom_number as eventtype -- distinct because only interested once in record of an event at a given time !!CHECK - i think wrong because there is the origin id which is unique per record
from 18_299_Lyratzopoulos_e2.cprd_clinical c
inner join 18_299_Lyratzopoulos.lookup_core_symptom l on c.medcode = l.medcode
inner join nadine.e2_allcrc l1 on l1.e_patid = c.e_patid
where l.symptom_number in (4, 14)
AND c.eventdate >= makedate(2011, 1) -- earliest can go back 1 year to 2011
AND c.eventdate >= l1.uts -- event within follow-up period
AND c.eventdate >= l1.crd_oneyear
AND  c.eventdate >= l1.age30_date
AND c.eventdate <= l1.lcd
AND (c.eventdate <= l1.tod or l1.tod is null)
AND (c.eventdate <= l1.deathdate or l1.deathdate is null)
AND c.eventdate < l1.age100_date
AND c.eventdate <= makedate(2018, 365)
);

select * from nadine.e2_allcrc_initial_sx_new limit 10;
select distinct year(eventdate) as sx_year from nadine.e2_allcrc_initial_sx_new order by sx_year; 
select min(eventdate), max(eventdate) from nadine.e2_allcrc_initial_sx_new; 
select count(*), count(distinct e_patid), eventtype from nadine.e2_allcrc_initial_sx_new group by eventtype;
select count(*), count(distinct e_patid) from nadine.e2_allcrc_initial_sx_new; 

-- Continue in R
-- Steps to be done:
	-- 3) Get eligible CRC
	-- 4) Keep symptoms for patients with eligilbe CRC only
	-- 5) Flag symptomatic CRC
	-- 6) Select only 1 CRC per patient (symptomatic or random)


select * from nadine.e2_allcrc_sx_final_dataset limit 10;
select count(*), count(distinct e_patid) from nadine.e2_allcrc_sx_final_dataset; 

select count(*), count(distinct e_patid) from nadine.e2_allcrc_sx_final_dataset 
where random_sample = 1 and cohort = 1 and stage_bin != "missing"; 

select count(*), count(distinct e_patid) from nadine.e2_allcrc_sx_final_dataset 
where random_sample = 1 and cohort = 0 and stage_bin != "missing"; 

select count(*), count(distinct e_patid) from nadine.e2_allcrc_sx_final_dataset 
where random_sample = 0 and cohort = 1 and stage_bin != "missing"; 

select count(*), count(distinct e_patid), symptomatic from nadine.e2_allcrc_sx_final_dataset 
where random_sample = 0 and cohort = 1 and stage_bin != "missing"
group by symptomatic; 

select count(*), count(distinct e_patid), symptomatic from nadine.e2_allcrc_sx_final_dataset 
where random_sample = 1 and cohort = 1 and stage_bin != "missing"
group by symptomatic;
