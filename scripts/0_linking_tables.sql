-- Purpose: Create first basic tables with all required information
-- Author:  Nadine Zakkak 

drop table if exists nadine.e2_cprd_case_file;
create table nadine.e2_cprd_case_file
select l1.e_pracid, p.*, l2.lcd, l2.uts, l3.imd2015_10
,case when l4.e_patid is not null then 1 else 0 end as cohort
,case when l5.e_patid is not null then 1 else 0 end as random_sample
from 18_299_Lyratzopoulos_e2.cprd_patient p
left join 18_299_Lyratzopoulos_e2.cprd_linkage_eligibility_gold l1 on p.e_patid = l1.e_patid
left join 18_299_Lyratzopoulos_e2.cprd_practice l2 on l1.e_pracid = l2.e_pracid
left join 18_299_Lyratzopoulos_e2.imd_2015 l3 on p.e_patid = l3.e_patid
left join 18_299_Lyratzopoulos_e2.cprd_cohort_file l4 on p.e_patid = l4.e_patid
left join 18_299_Lyratzopoulos_e2.cprd_random_sample l5 on p.e_patid = l5.e_patid;

select count(*), count(distinct e_patid) from 18_299_Lyratzopoulos_e2.cprd_patient; -- 2,871,070 both
select count(*), count(distinct e_patid), count(distinct e_pracid) from 18_299_Lyratzopoulos_e2.cprd_linkage_eligibility_gold; -- 2871070 patients and 408 practices
select count(*), count(distinct e_patid), count(distinct e_pracid) from nadine.e2_cprd_case_file; -- 2871070 patients and 408 practices
select count(*), count(distinct e_patid) from nadine.e2_cprd_case_file where random_sample = 1; -- 1m both
select count(*), count(distinct e_patid) from nadine.e2_cprd_case_file where cohort = 1; -- 2,530,253	both


alter table nadine.e2_cprd_case_file 
add column crd_oneyear date,
add column dob date,
add column age30_date date,
add column age100_date date;
set sql_safe_updates = 0;
update nadine.e2_cprd_case_file set crd_oneyear = date_add(crd, interval 1 year);
update nadine.e2_cprd_case_file set dob = makedate(yob, (365/2)); -- dob set to 'middle of the year'
update nadine.e2_cprd_case_file set age30_date = date_add(dob, interval 30 year);
update nadine.e2_cprd_case_file set age100_date = date_add(dob, interval 100 year);
set sql_safe_updates = 1;

CREATE INDEX `e_patid` ON nadine.e2_cprd_case_file (`e_patid`);
select * from nadine.e2_cprd_case_file limit 10;
select crd, crd_oneyear, yob, dob, age30_date, age100_date, deathdate, tod from nadine.e2_cprd_case_file limit 50;

select * from lookup_tables.lookup_cancersite;
select * from 18_299_Lyratzopoulos.lookup_core_cancersite;


