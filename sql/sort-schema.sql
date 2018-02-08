-- sort prototype schema

-- DROPS

drop view if exists perfect_variants_with_kits_assignments_and_unk;
drop view if exists perfect_variants_with_kits_assignments;
drop view if exists perfect_variants_with_kits;
drop view if exists perfect_variants;
drop view if exists kits_view;

drop table if exists s_variants;
drop table if exists s_calls;
-- drop table if exists s_call_passes;
-- drop table if exists s_call_fails;
-- drop table if exists s_sort_variants;
drop table if exists s_kits;
-- drop table if exists s_sort_kits;
drop table if exists s_dupes;
drop table if exists s_dupe_joins;

-- CREATES

create table s_variants (
 -- variant_id int, -- not needed for prototype
 variant_id int,  -- PK
 variant_loc varchar(10),
 name varchar(20) --,
 -- old_reference varchar(2), -- commenting out right now cuz not part of ian's doc
 -- sort_order int
);

create table s_calls(
 -- call_id int, -- PK - commenting out for now
 -- kit_id int, 
 kit_id varchar(10), 
 -- variant_loc int,
 variant_loc varchar(10),
 assigned boolean 
);

-- hide-me {{{

-- unique index (kit_id,variant_loc) -- commented out for now
/* (hide this for now){{{
-- note: I could see calls becoming maybe two tables.

create table call_passes (
 call_id int, -- PK
 variant_loc int -- PK
 -- new_reference varchar(2), -- not needed for prototype
 -- override int -- keep out for now 
);

create table call_fails (
 call_id int, -- PK
 variant_loc int --
 -- override int -- keep out for now   
);

}}} */

-- create table s_variants (
--  variant_loc int,
--  sort_order int
-- );

-- }}}

create table s_kits(
 -- kit_id  int,  -- later this can be person_id
 kit_id  varchar(10),  -- later this can be person_id
 sort_order int
);

-- hide-me {{{

-- create table s_sort_problems (
--  problem_id int -- PK
--  kit_id int, -- later can be person tbl
--  variants text,
--  issue_type tinyint,
--  notes text,
--  resolution text
-- );

-- create table s_call_overrides (
--   override_id -- PK  
--   call_id int -- 
--   old -- not sure how this could look
-- );

-- }}}

create table s_dupes (
  variant_id int, -- parent
  name varchar(30)
);

create table s_dupe_joins (
  variant_id int, -- parent
  dupe_variant_id int -- children
);

-- create table one_rec (foo int);
-- insert into one_rec (foo) values(1);

create view kits_view AS 
  SELECT DISTINCT kit_id from s_calls;

create view perfect_variants AS
  SELECT DISTINCT V.name, V.variant_loc, V.variant_id
  FROM s_calls C, s_variants V
  WHERE (C.assigned = -1 OR V.name = 'top') AND
  V.variant_loc = C.variant_loc;

create view perfect_variants_with_kits AS
  SELECT K.kit_id, PV.name, PV.variant_loc, PV.variant_id
  FROM perfect_variants PV
  CROSS JOIN 
  kits_view K;

create view perfect_variants_with_kits_assignments AS
  SELECT C.kit_id, PV.name, PV.variant_loc, PV.variant_id, C.assigned
  FROM s_calls C, perfect_variants PV
  WHERE C.variant_loc = PV.variant_loc;

create view perfect_variants_with_kits_assignments_and_unk AS
  SELECT PVK.kit_id, PVK.name, PVK.variant_loc, PVK.variant_id, ifnull(PVKA.assigned,0)
  FROM 
  perfect_variants_with_kits PVK
  LEFT JOIN 
  perfect_variants_with_kits_assignments PVKA
  ON
  PVK.variant_loc = PVKA.variant_loc AND
  PVK.kit_id = PVKA.kit_id;
