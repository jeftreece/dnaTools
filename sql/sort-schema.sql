-- sort prototype schema

-- DROP VIEWS {{{

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

drop view if exists get_max_snpnames;
drop view if exists x_perfect_variants_base;
drop view if exists x_perfect_variants_5;

-- drop view if exists x_perfect_assignments_with_unk_cnt_pos_v;
-- drop view if exists x_perfect_assignments_with_unk_cnt_neg_v;
-- drop view if exists x_perfect_assignments_with_unk_cnt_unk_v;

-- drop view if exists x_perfect_assignments_with_unk_cnt_pos_k;
-- drop view if exists x_perfect_assignments_with_unk_cnt_neg_k;
-- drop view if exists x_perfect_assignments_with_unk_cnt_unk_k;
 
drop view if exists x_pos_call_chk;
drop view if exists x_neg_call_chk;

drop view if exists get_max_snpnames;

-- }}}
-- DROP TABLES {{{

drop table if exists x_mx_kits;
drop table if exists x_mx_variants;
drop table if exists x_mx_idxs;
drop table if exists x_mx_calls;

-- }}}
-- CREATE TABLES {{{

create table x_mx_kits(
 ID int,
 kitId text
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
-- CREATE VIEWS (perfect) {{{

create view get_max_snpnames AS
 select max(snpname) as snpname,vID from snpnames group by 2;

create view x_kit_view AS 
  SELECT DISTINCT C.pID,max(D.kitId) as kitId from vcfcalls C, dataset D 
  WHERE C.pID=D.DNAID -- AND 
  -- C.pid IN (1049,1321,1707,1747,1801) 
  GROUP BY 1;

create view x_pos_call_chk as 
  SELECT DISTINCT vID from vcfcalls where assigned = 1;

create view x_neg_call_chk as 
  SELECT DISTINCT vID from vcfcalls where assigned = -1;

create view x_perfect_variants_base AS
  SELECT DISTINCT ifnull(S.snpname,V.ID) as name, V.pos, V.ID -- what is the right thing here? ID or pos with something?
  FROM vcfcalls C, x_pos_call_chk P, x_neg_call_chk N, variants V
  -- FROM vcfcalls C, variants V
  LEFT JOIN get_max_snpnames S
  ON S.vID = V.ID
  -- WHERE (C.assigned = -1 OR V.ID = -999) AND -- C.assigned
  WHERE -- C.assigned
  N.vID = V.ID and P.vID = V.ID AND 
  V.ID = C.vID and V.pos IN
  (13668461,7378685,12060401,19538924); -- z156, z381, z301, z28 -- u106, l48, z156, z8

create view x_perfect_variants_5 AS
  SELECT DISTINCT ifnull(S.snpname,V.ID) as name, V.pos, V.ID -- what is the right thing here? ID or pos with something?
  FROM vcfcalls C, x_pos_call_chk P, x_neg_call_chk N, variants V
  -- FROM vcfcalls C, variants V
  LEFT JOIN get_max_snpnames S
  ON S.vID = V.ID
  -- WHERE (C.assigned = -1 OR V.ID = -999) AND -- C.assigned
  WHERE -- C.assigned
  N.vID = V.ID and P.vID = V.ID AND 
  V.ID = C.vID and 
  V.pos not in (13668461,7378685,12060401,19538924) 
  LIMIT 5;

create view x_perfect_variants AS
  SELECT * from x_perfect_variants_base
  UNION
  SELECT * from x_perfect_variants_5;

create view x_perfect_variants_with_kits AS
  SELECT K.pID, PV.name, PV.pos, PV.ID as vID,K.kitId
  FROM x_perfect_variants PV
  CROSS JOIN x_kit_view K;

create view x_perfect_assignments AS
  SELECT C.pID, PV.name, PV.pos, PV.ID as vID, C.assigned -- C.assigned
  FROM vcfcalls C, x_perfect_variants PV
  WHERE C.vID = PV.ID;

create view x_perfect_assignments_with_unk AS
  SELECT PVK.pID, PVK.name, ifnull(PVKA.assigned,0) as assigned, PVK.pos, PVK.vID,PVK.kitId
  FROM x_perfect_variants_with_kits PVK 
  LEFT JOIN x_perfect_assignments PVKA
  ON PVK.vID = PVKA.vID AND
  PVK.pID = PVKA.pID;

-- create view x_perfect_assignments_with_unk_cnt_pos_v AS
--  SELECT sum(case when assigned <> 1 then 0 else 1 end) as cnt,name,assigned,vID from x_perfect_assignments_with_unk 
--  GROUP BY 4;

-- create view x_perfect_assignments_with_unk_cnt_neg_v AS
--  SELECT sum(case when assigned <> -1 then 0 else 1 end) as cnt,name,assigned,vID from x_perfect_assignments_with_unk 
--  GROUP BY 4;

-- create view x_perfect_assignments_with_unk_cnt_unk_v AS
--  SELECT sum(case when assigned <> 0 then 0 else 1 end) as cnt,name,assigned,vID from x_perfect_assignments_with_unk 
--  GROUP BY 4;

-- create view x_perfect_assignments_with_unk_cnt_pos_k AS
--  SELECT sum(case when assigned <> 1 then 0 else 1 end) as cnt,pID,assigned from x_perfect_assignments_with_unk 
--  GROUP BY 2;

-- create view x_perfect_assignments_with_unk_cnt_neg_k AS
--  SELECT sum(case when assigned <> -1 then 0 else 1 end) as cnt,pID,assigned from x_perfect_assignments_with_unk 
--  GROUP BY 2;

-- create view x_perfect_assignments_with_unk_cnt_unk_k AS
--  SELECT sum(case when assigned <> 0 then 0 else 1 end) as cnt,pID,assigned from x_perfect_assignments_with_unk 
--  GROUP BY 2;

-- }}}
-- CREATE VIEWS (saved) {{{

create view x_saved_kits AS 
  SELECT DISTINCT ID,kitId from x_mx_kits;

create view x_saved_variants AS 
  SELECT DISTINCT name, pos, ID
  FROM x_mx_variants;

create view x_saved_variants_with_kits AS
  SELECT SK.ID as pID, SV.name, SV.pos, SV.ID as vID,SK.kitID
  FROM x_saved_variants SV
  CROSS JOIN x_saved_kits SK;

create view x_saved_assignments AS
  SELECT SC.pID, SV.name, SV.pos, SV.ID as vID, SC.assigned
  FROM x_mx_calls SC, x_saved_variants SV
  WHERE SC.vID = SV.ID;

create view x_saved_assignments_with_unk AS
  SELECT SVK.pID, SVK.name, ifnull(SVKA.assigned,0) as assigned, SVK.pos, SVK.vID, SVK.kitId, VI.idx, KI.idx
  FROM x_mx_idxs VI, x_mx_idxs KI, x_saved_variants_with_kits SVK
  LEFT JOIN x_saved_assignments SVKA
  ON SVK.vID = SVKA.vID AND SVK.pID = SVKA.pID
  -- WHERE VI.type_id = 0 AND VI.axis_id = SVK.vID AND  -- fix this
  WHERE VI.type_id = 0 AND VI.axis_id = SVK.name AND  -- fix this
  KI.type_id = 1 AND KI.axis_id = SVK.kitId
  ORDER BY 7,8;

-- }}}

