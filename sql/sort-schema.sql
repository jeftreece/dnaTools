-- sort prototype schema

-- DROP VIEWS

drop view if exists saved_assignments_with_unk;
drop view if exists saved_assignments;
drop view if exists saved_variants_with_kits;
drop view if exists saved_variants;
drop view if exists saved_kits;

drop view if exists perfect_assignments_with_unk;
drop view if exists perfect_assignments;
drop view if exists perfect_variants_with_kits;
drop view if exists perfect_variants;
drop view if exists kits_view;

-- DROP TABLES

drop table if exists s_variants;
drop table if exists s_calls;
-- drop table if exists s_call_passes;
-- drop table if exists s_call_fails;
-- drop table if exists s_sort_variants;
drop table if exists s_kits;
-- drop table if exists s_sort_kits;
drop table if exists s_dupes;
drop table if exists s_dupe_joins;

drop table if exists s_mx_idxs;
drop table if exists s_mx_variants;
drop table if exists s_mx_calls;

-- CREATE TABLES

create table s_variants (
 -- variant_id int, -- not needed for prototype
 variant_id int,  -- PK
 variant_loc varchar(10),
 name varchar(20) --,
 -- old_reference varchar(2), -- commenting out right now cuz not part of ian's doc
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
 kit_id  varchar(10)  -- later this can be person_id
);

create table s_mx_kits(
 -- kit_id  int,  -- later this can be person_id
 kit_id  varchar(10),  -- later this can be person_id
 idx int
);

create table s_mx_variants (
 -- variant_id int, -- not needed for prototype
 variant_id int,  -- PK
 ref_variant_id int,
 variant_loc varchar(10),
 name varchar(20), --,
 idx int
 -- old_reference varchar(2), -- commenting out right now cuz not part of ian's doc
);

-- create table s_mx_idxs(
-- type_id int,           -- 0 = variants, 1 = kits
-- axis_id varchar(10),   -- either the variant id or the kit_id 
--);

create table s_mx_calls (
 kit_id varchar(10),        
 variant_loc varchar(10),
 assigned boolean 
 confidence int,
 changed int
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

-- CREATE VIEWS

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
  CROSS JOIN kits_view K;

create view perfect_assignments AS
  SELECT C.kit_id, PV.name, PV.variant_loc, PV.variant_id, C.assigned
  FROM s_calls C, perfect_variants PV
  WHERE C.variant_loc = PV.variant_loc;

create view perfect_assignments_with_unk AS
  SELECT PVK.kit_id, PVK.name, ifnull(PVKA.assigned,0), PVK.variant_loc, PVK.variant_id
  FROM perfect_variants_with_kits PVK
  LEFT JOIN perfect_assignments PVKA
  ON PVK.variant_loc = PVKA.variant_loc AND
  PVK.kit_id = PVKA.kit_id;

create view saved_kits AS 
  SELECT DISTINCT kit_id, idx from s_mx_kits;

create view saved_variants AS 
  SELECT DISTINCT name, variant_loc, variant_id, idx
  FROM s_mx_variants;

create view saved_variants_with_kits AS
  SELECT SK.kit_id, SV.name, SV.variant_loc, SV.variant_id
  FROM saved_variants SV
  CROSS JOIN saved_kits SK;

create view saved_assignments AS
  SELECT SC.kit_id, SV.name, SV.variant_loc, SV.variant_id, SC.assigned
  FROM s_mx_calls SC, saved_variants SV
  WHERE SC.variant_loc = SV.variant_loc;

create view saved_assignments_with_unk AS
  SELECT SVK.kit_id, SVK.name, ifnull(SVKA.assigned,0), SVK.variant_loc, SVK.variant_id
  FROM saved_variants_with_kits SVK
  LEFT JOIN saved_assignments SVKA
  ON SVK.variant_loc = SVKA.variant_loc AND
  SVK.kit_id = SVKA.kit_id;

