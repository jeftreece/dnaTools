-- sort prototype schema

-- DROP VIEWS {{{

DROP VIEW IF EXISTS v_imx_kits;
DROP VIEW IF EXISTS v_imx_variants;
DROP VIEW IF EXISTS v_imx_variants_test;
DROP VIEW IF EXISTS v_imx_variants_all;
DROP VIEW IF EXISTS v_imx_variants_with_kits;
DROP VIEW IF EXISTS v_imx_assignments;
DROP VIEW IF EXISTS v_imx_assignments_with_unk;

DROP VIEW IF EXISTS v_mx_kits;
DROP VIEW IF EXISTS v_mx_variants;

-- DROP VIEW IF EXISTS x_saved_variants_with_kits;
-- DROP VIEW IF EXISTS x_saved_assignments;
-- DROP VIEW IF EXISTS x_saved_assignments_with_unk;

DROP VIEW IF EXISTS v_imx_variants_base;
DROP VIEW IF EXISTS v_imx_variants_lim;

DROP VIEW IF EXISTS v_pos_call_chk;
DROP VIEW IF EXISTS v_neg_call_chk;
DROP VIEW IF EXISTS v_max_snpnames;

DROP VIEW IF EXISTS v_only_pos_variants;
DROP VIEW IF EXISTS v_only_neg_variants;

-- }}}
-- DROP TABLES {{{

DROP TABLE IF EXISTS mx_kits;
DROP TABLE IF EXISTS mx_variants;
DROP TABLE IF EXISTS mx_idxs;
DROP TABLE IF EXISTS mx_calls;
DROP TABLE IF EXISTS mx_dupe_variants;

-- }}}
-- DROP INDEXES {{{

DROP INDEX IF EXISTS snpidx;

-- }}}
-- CREATE TABLES {{{

CREATE TABLE mx_kits(
 ID int,
 kitId text
);

CREATE TABLE mx_variants (
 ID int,  
 ref_variant_id int,
 name text,
 pos int
);

CREATE TABLE mx_dupe_variants (
 vId int,  
 dupe_vID int
);

CREATE TABLE mx_idxs(
 type_id int,   -- 0 = variants, 1 = kits (people)
 axis_id int,   -- either the variant id (vID) or the kit_id (pID) 
 idx int
);

CREATE TABLE mx_calls (
 pID int,
 vID int,
 assigned boolean,
 confidence int,
 changed int
);

-- }}}
-- CREATE VIEWS (import to matrix) {{{

CREATE VIEW v_max_snpnames AS
  SELECT max(snpname) as snpname,vID from snpnames group by 2;

CREATE VIEW v_imx_kits AS 
  SELECT DISTINCT C.pID,max(D.kitId) as kitId from vcfcalls C, dataset D 
  WHERE C.pID=D.ID 
  GROUP BY 1;

CREATE VIEW v_pos_call_chk as 
  SELECT DISTINCT vID from vcfcalls where assigned = 1;

CREATE VIEW v_neg_call_chk as 
  SELECT DISTINCT vID from vcfcalls where assigned = -1;

-- To create test data
-- --------------------
CREATE VIEW v_imx_variants_base AS
  SELECT DISTINCT ifnull(S.snpname,V.ID) as name, V.pos, V.ID 
  FROM vcfcalls C, v_pos_call_chk P, v_neg_call_chk N, variants V 
  LEFT JOIN v_max_snpnames S
  ON S.vID = V.ID
  WHERE
  N.vID = V.ID AND P.vID = V.ID AND 
  V.ID = C.vID AND V.pos IN
  (3019783,15732138,20577481,8928037,21450311,6920349,12879820,13668461,19995425,20029258,7378686,12060401,19538924,20323911);

-- To create test data
-- --------------------
CREATE VIEW v_imx_variants_lim AS
  SELECT DISTINCT ifnull(S.snpname,V.ID) as name, V.pos, V.ID 
  FROM vcfcalls C, v_pos_call_chk P, v_neg_call_chk N, variants V
  LEFT JOIN v_max_snpnames S
  ON S.vID = V.ID
  WHERE 
  N.vID = V.ID AND P.vID = V.ID AND
  V.ID = C.vID AND 
  V.pos NOT IN 
  (3019783,15732138,20577481,8928037,21450311,6920349,12879820,13668461,19995425,20029258,7378686,12060401,19538924,20323911)
  LIMIT 20;

-- To create test data
-- --------------------
CREATE VIEW v_imx_variants_test AS
  SELECT * from v_imx_variants_base
  UNION
  SELECT * from v_imx_variants_lim;

-- To create all data
-- --------------------
CREATE VIEW v_imx_variants_all AS
  SELECT DISTINCT ifnull(S.snpname,V.ID) as name, V.pos, V.ID 
  FROM vcfcalls C, v_pos_call_chk P, v_neg_call_chk N, variants V
  LEFT JOIN v_max_snpnames S
  ON S.vID = V.ID
  WHERE 
  N.vID = V.ID AND P.vID = V.ID AND
  V.ID = C.vID;

CREATE VIEW v_imx_variants AS
  SELECT * from v_imx_variants_all; -- this is a toggle btw test and all

CREATE VIEW v_only_pos_variants AS
  SELECT DISTINCT P.vID from v_pos_call_chk P 
  LEFT JOIN v_neg_call_chk N
  ON P.vID = N.vID
  WHERE N.vID is Null;
  
CREATE VIEW v_only_neg_variants AS
  SELECT DISTINCT N.vID from v_neg_call_chk N
  LEFT JOIN v_pos_call_chk P
  ON P.vID = N.vID
  WHERE P.vID is Null;
  
CREATE VIEW v_imx_variants_with_kits AS
  SELECT DISTINCT K.pID, PV.name, PV.pos, PV.ID as vID, K.kitId
  FROM v_imx_variants PV
  CROSS JOIN v_imx_kits K;

CREATE VIEW v_imx_assignments AS
  SELECT DISTINCT C.pID, PV.name, PV.pos, PV.ID as vID, C.assigned 
  FROM vcfcalls C, v_imx_variants PV
  WHERE C.vID = PV.ID;

CREATE VIEW v_imx_assignments_with_unk AS
  SELECT DISTINCT PVK.pID, PVK.name, ifnull(PVKA.assigned,0) as assigned, PVK.pos, PVK.vID, PVK.kitId
  FROM v_imx_variants_with_kits PVK 
  LEFT JOIN v_imx_assignments PVKA
  ON PVK.vID = PVKA.vID AND
  PVK.pID = PVKA.pID;

-- }}}
-- CREATE VIEWS (inside matrix) {{{

CREATE VIEW v_mx_kits AS 
  SELECT DISTINCT ID,kitId from mx_kits;

CREATE VIEW v_mx_variants AS 
  SELECT DISTINCT name, pos, ID
  FROM mx_variants;

-- }}}
-- CREATE INDEXES{{{

CREATE INDEX snpidx on snpnames(snpname);

/*}}}*/
