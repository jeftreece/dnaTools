-- sort prototype schema

-- (old) DROP VIEWS{{{

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

-- }}}
-- (old) DROP TABLES {{{

 drop table if exists s_variants;
 drop table if exists s_calls;
 drop table if exists s_kits;

 drop table if exists s_mx_kits;
 drop table if exists s_mx_idxs;
 drop table if exists s_mx_variants;
 drop table if exists s_mx_calls;

-- }}}
-- (old) CREATE TABLES {{{

-- create table s_variants (
--  variant_id int,  
--  variant_loc varchar(10),
--  name varchar(20)
-- );

-- create table s_calls(
--  kit_id varchar(10), 
--  variant_loc varchar(10),
--  assigned boolean 
-- );

-- create table s_kits(
--  kit_id  varchar(10)
-- );

-- create table s_mx_kits(
--  kit_id  varchar(10),
--  idx int
-- );

-- create table s_mx_variants (
--  variant_id int,  
--  ref_variant_id int,
--  variant_loc varchar(10),
--  name varchar(20), 
--  idx int
-- );

-- create table s_mx_idxs(
--  type_id int,           -- 0 = variants, 1 = kits
--  axis_id varchar(10),   -- either the variant id or the kit_id 
--  idx int
-- );

-- create table s_mx_calls (
--  kit_id varchar(10),        
--  variant_loc varchar(10),
--  assigned boolean,
--  confidence int,
--  changed int
-- );

-- }}}
-- (old) CREATE VIEWS (perfect) {{{

-- create view kits_view AS 
--   SELECT DISTINCT kit_id from s_calls;

-- create view perfect_variants AS
--   SELECT DISTINCT V.name, V.variant_loc, V.variant_id
--   FROM s_calls C, s_variants V
--   WHERE (C.assigned = -1 OR V.name = 'top') AND
--   V.variant_loc = C.variant_loc;

-- create view perfect_variants_with_kits AS
--   SELECT K.kit_id, PV.name, PV.variant_loc, PV.variant_id
--   FROM perfect_variants PV
--   CROSS JOIN kits_view K;

-- create view perfect_assignments AS
--   SELECT C.kit_id, PV.name, PV.variant_loc, PV.variant_id, C.assigned
--   FROM s_calls C, perfect_variants PV
--   WHERE C.variant_loc = PV.variant_loc;

-- create view perfect_assignments_with_unk AS
--   SELECT PVK.kit_id, PVK.name, ifnull(PVKA.assigned,0) as assigned, PVK.variant_loc, PVK.variant_id
--   FROM perfect_variants_with_kits PVK
--   LEFT JOIN perfect_assignments PVKA
--   ON PVK.variant_loc = PVKA.variant_loc AND
--   PVK.kit_id = PVKA.kit_id;

-- }}}
-- (old) CREATE VIEWS (saved) {{{

-- create view saved_kits AS 
--   SELECT DISTINCT kit_id from s_mx_kits;

-- create view saved_variants AS 
--   SELECT DISTINCT name, variant_loc, variant_id
--   FROM s_mx_variants;

-- create view saved_variants_with_kits AS
--   SELECT SK.kit_id, SV.name, SV.variant_loc, SV.variant_id
--   FROM saved_variants SV
--   CROSS JOIN saved_kits SK;

-- create view saved_assignments AS
--   SELECT SC.kit_id, SV.name, SV.variant_loc, SV.variant_id, SC.assigned
--   FROM s_mx_calls SC, saved_variants SV
--   WHERE SC.variant_loc = SV.variant_loc;

-- create view saved_assignments_with_unk AS
--   SELECT SVK.kit_id, SVK.name, ifnull(SVKA.assigned,0) as assigned, SVK.variant_loc, SVK.variant_id, VI.idx, KI.idx
--   FROM s_mx_idxs VI, s_mx_idxs KI, saved_variants_with_kits SVK
--   LEFT JOIN saved_assignments SVKA
--   ON SVK.variant_loc = SVKA.variant_loc AND
--   SVK.kit_id = SVKA.kit_id
--   WHERE VI.type_id = 0 and VI.axis_id = SVK.name AND 
--   KI.type_id = 1 and KI.axis_id = SVK.kit_id
--   ORDER BY 6,7;

-- }}}

-- (new) DROP VIEWS {{{

drop view if exists x_kit_view;
drop view if exists x_perfect_variants;
drop view if exists x_perfect_variants_with_kits;
drop view if exists x_perfect_assignments;
drop view if exists x_perfect_assignments_with_unk;

drop view if exists x_saved_kits;
drop view if exists x_saved_variants;
drop view if exists x_saved_variants_with_kits;
drop view if exists x_saved_assignments;
drop view if exists x_saved_assignments_with_unk;

drop view if exists x_perfect_assignments_with_unk_cnt_pos_v;
drop view if exists x_perfect_assignments_with_unk_cnt_neg_v;
drop view if exists x_perfect_assignments_with_unk_cnt_unk_v;

drop view if exists x_perfect_assignments_with_unk_cnt_pos_k;
drop view if exists x_perfect_assignments_with_unk_cnt_neg_k;
drop view if exists x_perfect_assignments_with_unk_cnt_unk_k;

-- }}}
-- (new) DROP TABLES {{{

drop table if exists x_mx_kits;
drop table if exists x_mx_variants;
drop table if exists x_mx_idxs;
drop table if exists x_mx_calls;

-- }}}
-- (new) CREATE TABLES {{{

create table x_mx_kits(
 ID int,
 idx int
);

create table x_mx_variants (
 ID int,  
 ref_variant_id int,
 pos int,
 name text
);

create table x_mx_idxs(
 type_id int,           -- 0 = variants, 1 = kits (people)
 axis_id int,   -- either the variant id (vID) or the kit_id (pID) 
 idx int
);

create table x_mx_calls (
 pID int,
 vID int,
 assigned boolean,
 confidence int,
 changed int
);

-- }}}
-- (new) CREATE VIEWS (perfect) {{{

create view x_kit_view AS 
  SELECT DISTINCT pID from vcfcalls;

-- create view x_perfect_variants AS
--   SELECT DISTINCT ifnull(S.snpname,V.pos) as name, V.pos, V.ID
--   FROM vcfcalls C, variants V 
--   LEFT JOIN snpnames S
--   ON S.vID = V.ID
--   WHERE (C.assigned = -1 OR V.ID = -999) AND -- C.assigned
--   V.ID = C.vID;

create view x_perfect_variants AS
  SELECT DISTINCT ifnull(S.snpname,V.ID) as name, V.pos, V.ID -- what is the right thing here? ID or pos with something?
  FROM vcfcalls C, variants V 
  LEFT JOIN snpnames S
  ON S.vID = V.ID
  WHERE (C.assigned = -1 OR V.ID = -999) AND -- C.assigned
  V.ID = C.vID limit 5;

create view x_perfect_variants_with_kits AS
  SELECT K.pID, PV.name, PV.pos, PV.ID as vID
  FROM x_perfect_variants PV
  CROSS JOIN x_kit_view K;

create view x_perfect_assignments AS
  SELECT C.pID, PV.name, PV.pos, PV.ID as vID, C.assigned -- C.assigned
  FROM vcfcalls C, x_perfect_variants PV
  WHERE C.vID = PV.ID;

create view x_perfect_assignments_with_unk AS
  SELECT PVK.pID, PVK.name, ifnull(PVKA.assigned,0) as assigned, PVK.pos, PVK.vID
  FROM x_perfect_variants_with_kits PVK 
  LEFT JOIN x_perfect_assignments PVKA
  ON PVK.vID = PVKA.vID AND
  PVK.pID = PVKA.pID;

create view x_perfect_assignments_with_unk_cnt_pos_v AS
 SELECT sum(case when assigned <> 1 then 0 else 1 end) as cnt,name,assigned,vID from x_perfect_assignments_with_unk 
 GROUP BY 4;

create view x_perfect_assignments_with_unk_cnt_neg_v AS
 SELECT sum(case when assigned <> -1 then 0 else 1 end) as cnt,name,assigned,vID from x_perfect_assignments_with_unk 
 GROUP BY 4;

create view x_perfect_assignments_with_unk_cnt_unk_v AS
 SELECT sum(case when assigned <> 0 then 0 else 1 end) as cnt,name,assigned,vID from x_perfect_assignments_with_unk 
 GROUP BY 4;

create view x_perfect_assignments_with_unk_cnt_pos_k AS
 SELECT sum(case when assigned <> 1 then 0 else 1 end) as cnt,pID,assigned from x_perfect_assignments_with_unk 
 GROUP BY 2;

create view x_perfect_assignments_with_unk_cnt_neg_k AS
 SELECT sum(case when assigned <> -1 then 0 else 1 end) as cnt,pID,assigned from x_perfect_assignments_with_unk 
 GROUP BY 2;

create view x_perfect_assignments_with_unk_cnt_unk_k AS
 SELECT sum(case when assigned <> 0 then 0 else 1 end) as cnt,pID,assigned from x_perfect_assignments_with_unk 
 GROUP BY 2;

-- }}}
-- (new) CREATE VIEWS (saved) {{{

create view x_saved_kits AS 
  SELECT DISTINCT ID from x_mx_kits;

create view x_saved_variants AS 
  SELECT DISTINCT name, pos, ID
  FROM x_mx_variants limit 5;

create view x_saved_variants_with_kits AS
  SELECT SK.ID as pID, SV.name, SV.pos, SV.ID as vID
  FROM x_saved_variants SV
  CROSS JOIN x_saved_kits SK;

create view x_saved_assignments AS
  SELECT SC.pID, SV.name, SV.pos, SV.ID as vID, SC.assigned
  FROM x_mx_calls SC, x_saved_variants SV
  WHERE SC.vID = SV.ID;

create view x_saved_assignments_with_unk AS
  SELECT SVK.pID, SVK.name, ifnull(SVKA.assigned,0) as assigned, SVK.pos, SVK.vID, VI.idx, KI.idx
  FROM x_mx_idxs VI, x_mx_idxs KI, x_saved_variants_with_kits SVK
  LEFT JOIN x_saved_assignments SVKA
  ON SVK.vID = SVKA.vID AND SVK.pID = SVKA.pID
  -- WHERE VI.type_id = 0 AND VI.axis_id = SVK.vID AND  -- fix this
  WHERE VI.type_id = 0 AND VI.axis_id = SVK.name AND  -- fix this
  KI.type_id = 1 AND KI.axis_id = SVK.pID
  ORDER BY 6,7;

-- }}}

