//genesis
// kkit Version 11 flat dumpfile
 
// Saved on Thu Sep 19 13:57:35 2024
 
include kkit {argv 1}
 
FASTDT = 0.001
SIMDT = 0.001
CONTROLDT = 0.01
PLOTDT = 0.001
MAXTIME = 100
TRANSIENT_TIME = 2
VARIABLE_DT_FLAG = 0
DEFAULT_VOL = 1.1587e-15
VERSION = 11.0
setfield /file/modpath value ~/scripts/modules
kparms
 
//genesis

initdump -version 3 -ignoreorphans 1
simobjdump doqcsinfo filename accessname accesstype transcriber developer \
  citation species tissue cellcompartment methodology sources \
  model_implementation model_validation x y z
simobjdump table input output alloced step_mode stepsize x y z
simobjdump xtree path script namemode sizescale
simobjdump xcoredraw xmin xmax ymin ymax
simobjdump xtext editable
simobjdump xgraph xmin xmax ymin ymax overlay
simobjdump xplot pixflags script fg ysquish do_slope wy
simobjdump group xtree_fg_req xtree_textfg_req plotfield expanded movealone \
  link savename file version md5sum mod_save_flag x y z
simobjdump geometry size dim shape outside xtree_fg_req xtree_textfg_req x y \
  z
simobjdump kpool DiffConst CoInit Co n nInit mwt nMin vol slave_enable \
  geomname xtree_fg_req xtree_textfg_req x y z
simobjdump kreac kf kb notes xtree_fg_req xtree_textfg_req x y z
simobjdump kenz CoComplexInit CoComplex nComplexInit nComplex vol k1 k2 k3 \
  keepconc usecomplex notes xtree_fg_req xtree_textfg_req link x y z
simobjdump stim level1 width1 delay1 level2 width2 delay2 baselevel trig_time \
  trig_mode notes xtree_fg_req xtree_textfg_req is_running x y z
simobjdump xtab input output alloced step_mode stepsize notes editfunc \
  xtree_fg_req xtree_textfg_req baselevel last_x last_y is_running x y z
simobjdump kchan perm gmax Vm is_active use_nernst notes xtree_fg_req \
  xtree_textfg_req x y z
simobjdump transport input output alloced step_mode stepsize dt delay clock \
  kf xtree_fg_req xtree_textfg_req x y z
simobjdump proto x y z
simobjdump text str
simundump geometry /kinetics/geometry 0 1.15874090871e-18 3 sphere "" white \
  black 4 2 0
simundump geometry /kinetics/geometry[1] 0 1.1587e-15 3 sphere "" white black \
  0 0 0
simundump text /kinetics/notes 0 ""
call /kinetics/notes LOAD \
""
simundump text /kinetics/geometry/notes 0 ""
call /kinetics/geometry/notes LOAD \
""
simundump text /kinetics/geometry[1]/notes 0 ""
call /kinetics/geometry[1]/notes LOAD \
""
simundump group /kinetics/CaM 0 blue green x 0 0 "" defaultfile defaultfile.g \
  0 0 0 1 2 0
simundump text /kinetics/CaM/notes 0 ""
call /kinetics/CaM/notes LOAD \
""
simundump kpool /kinetics/CaM/CaM_Ca4 1 2e-11 0 0 0 0 0 0 600 0 \
  /kinetics/geometry blue 10 2551 5857 0
simundump text /kinetics/CaM/CaM_Ca4/notes 0 ""
call /kinetics/CaM/CaM_Ca4/notes LOAD \
""
simundump kpool /kinetics/CaM/CaM 1 2e-11 2 2 1200 1200 0 0 600 0 \
  /kinetics/geometry pink 10 2182 5857 0
simundump text /kinetics/CaM/CaM/notes 0 \
  "There is a LOT of this in the cell: upto 1% of total protein mass. (Alberts et al) Say 25 uM. Meyer et al Science 256 1199-1202 1992 refer to studies saying it is comparable to CaMK levels."
call /kinetics/CaM/CaM/notes LOAD \
"There is a LOT of this in the cell: upto 1% of total protein mass. (Alberts et al) Say 25 uM. Meyer et al Science 256 1199-1202 1992 refer to studies saying it is comparable to CaMK levels."
simundump kpool /kinetics/CaM/CaM_Ca3 1 2e-11 0 0 0 0 0 0 600 0 \
  /kinetics/geometry hotpink 10 3280 5857 0
simundump text /kinetics/CaM/CaM_Ca3/notes 0 ""
call /kinetics/CaM/CaM_Ca3/notes LOAD \
""
simundump kpool /kinetics/CaM/CaM_Ca2 1 2e-11 0 0 0 0 0 0 600 0 \
  /kinetics/geometry pink 10 2471 7014 0
simundump text /kinetics/CaM/CaM_Ca2/notes 0 \
  "This is the intermediate where the TR2 end (the high-affinity end) has bound the Ca but the TR1 end has not."
call /kinetics/CaM/CaM_Ca2/notes LOAD \
"This is the intermediate where the TR2 end (the high-affinity end) has bound the Ca but the TR1 end has not."
simundump kpool /kinetics/CaM/CaM_Ca 1 2e-11 0 0 0 0 0 0 600 0 \
  /kinetics/geometry pink 10 3201 7014 0
simundump text /kinetics/CaM/CaM_Ca/notes 0 \
  "This is the intermediate where the TR2 end (the high-affinity end) has bound the Ca but the TR1 end has not."
call /kinetics/CaM/CaM_Ca/notes LOAD \
"This is the intermediate where the TR2 end (the high-affinity end) has bound the Ca but the TR1 end has not."
simundump kreac /kinetics/CaM/CaM_Ca3_bind_Ca 0 0.003 10 "" white 10 3004 \
  6478 0
simundump text /kinetics/CaM/CaM_Ca3_bind_Ca/notes 0 \
  "Use K3 = 21.5 uM here from Stemmer and Klee table 3. kb/kf =21.5 * 6e5 so kf = 7.75e-7, kb = 10"
call /kinetics/CaM/CaM_Ca3_bind_Ca/notes LOAD \
"Use K3 = 21.5 uM here from Stemmer and Klee table 3. kb/kf =21.5 * 6e5 so kf = 7.75e-7, kb = 10"
simundump kreac /kinetics/CaM/CaM_bind_Ca 0 0.014141 8.4853 "" white 10 3333 \
  6473 0
simundump text /kinetics/CaM/CaM_bind_Ca/notes 0 \
  "Lets use the fast rate consts here. Since the rates are so different, I am not sure whether the order is relevant. These correspond to the TR2C fragment. We use the Martin et al rates here, plus the Drabicowski binding consts. All are scaled by 3X to cell temp. kf = 2e-10 kb = 72 Stemmer & Klee: K1=.9, K2=1.1. Assume 1.0uM for both. kb/kf=3.6e11. If kb=72, kf = 2e-10 (Exactly the same !)...."
call /kinetics/CaM/CaM_bind_Ca/notes LOAD \
"Lets use the fast rate consts here. Since the rates are so different, I am not sure whether the order is relevant. These correspond to the TR2C fragment. We use the Martin et al rates here, plus the Drabicowski binding consts. All are scaled by 3X to cell temp. kf = 2e-10 kb = 72 Stemmer & Klee: K1=.9, K2=1.1. Assume 1.0uM for both. kb/kf=3.6e11. If kb=72, kf = 2e-10 (Exactly the same !) 19 May 2006. Splitting the old CaM-TR2-bind-Ca reaction into two steps, each binding 1 Ca. This improves numerical stability and is conceptually better too. Overall rates are the same, so each kf and kb is the square root of the earlier ones. So kf = 1.125e-4, kb = 8.4853"
simundump kreac /kinetics/CaM/CaM_Ca2_bind_Ca 0 0.006 10 "" white 10 3667 \
  6473 0
simundump text /kinetics/CaM/CaM_Ca2_bind_Ca/notes 0 \
  "K3 = 21.5, K4 = 2.8. Assuming that the K4 step happens first, we get kb/kf = 2.8 uM = 1.68e6 so kf =6e-6 assuming kb = 10"
call /kinetics/CaM/CaM_Ca2_bind_Ca/notes LOAD \
"K3 = 21.5, K4 = 2.8. Assuming that the K4 step happens first, we get kb/kf = 2.8 uM = 1.68e6 so kf =6e-6 assuming kb = 10"
simundump kreac /kinetics/CaM/CaM_Ca_bind_Ca 0 0.014141 8.4853 "" white 10 \
  4000 6455 0
simundump text /kinetics/CaM/CaM_Ca_bind_Ca/notes 0 \
  "Lets use the fast rate consts here. Since the rates are so different, I am not sure whether the order is relevant. These correspond to the TR2C fragment. We use the Martin et al rates here, plus the Drabicowski binding consts. All are scaled by 3X to cell temp. kf = 2e-10 kb = 72 Stemmer & Klee: K1=.9, K2=1.1. Assume 1.0uM for both. kb/kf=3.6e11. If kb=72, kf = 2e-10 (Exactly the same !)...."
call /kinetics/CaM/CaM_Ca_bind_Ca/notes LOAD \
"Lets use the fast rate consts here. Since the rates are so different, I am not sure whether the order is relevant. These correspond to the TR2C fragment. We use the Martin et al rates here, plus the Drabicowski binding consts. All are scaled by 3X to cell temp. kf = 2e-10 kb = 72 Stemmer & Klee: K1=.9, K2=1.1. Assume 1.0uM for both. kb/kf=3.6e11. If kb=72, kf = 2e-10 (Exactly the same !) 19 May 2006. Splitting the old CaM-TR2-bind-Ca reaction into two steps, each binding 1 Ca. This improves numerical stability and is conceptually better too. Overall rates are the same, so each kf and kb is the square root of the earlier ones. So kf = 1.125e-4, kb = 8.4853"
simundump group /kinetics/Ca 0 blue green x 0 0 "" defaultfile defaultfile.g \
  0 0 0 1 2 0
simundump text /kinetics/Ca/notes 0 ""
call /kinetics/Ca/notes LOAD \
""
simundump kpool /kinetics/Ca/Ca 0 1e-10 0.08 0.08 48 48 0 0 600 0 \
  /kinetics/geometry red 18 3513 7522 0
simundump text /kinetics/Ca/Ca/notes 0 ""
call /kinetics/Ca/Ca/notes LOAD \
""
simundump kpool /kinetics/Ca/Ca_input 0 0 0 0 0 0 0 0 600 4 \
  /kinetics/geometry orange 18 2513 7522 0
simundump text /kinetics/Ca/Ca_input/notes 0 ""
call /kinetics/Ca/Ca_input/notes LOAD \
""
simundump kreac /kinetics/Ca/Ca_stim 0 20 0 "" white 40 3013 7522 0
simundump text /kinetics/Ca/Ca_stim/notes 0 ""
call /kinetics/Ca/Ca_stim/notes LOAD \
""
simundump kpool /kinetics/Ca/Ca_ext 0 0 0.08 0.08 48 48 0 0 600 4 \
  /kinetics/geometry[1] blue black 3229 7921 0
simundump text /kinetics/Ca/Ca_ext/notes 0 ""
call /kinetics/Ca/Ca_ext/notes LOAD \
""
simundump kreac /kinetics/Ca/Ca_homeostasis 0 0.5 0.5 "" white black 3317 \
  7644 0
simundump text /kinetics/Ca/Ca_homeostasis/notes 0 ""
call /kinetics/Ca/Ca_homeostasis/notes LOAD \
""
simundump xgraph /graphs/conc1 0 0 99 0.001 0.999 0
simundump xgraph /graphs/conc2 0 0 100 0 1 0
simundump xplot /graphs/conc1/Ca.Co 3 524288 \
  "delete_plot.w <s> <d>; edit_plot.D <w>" red 0 0 1
simundump xplot /graphs/conc1/CaM_Ca4.Co 3 524288 \
  "delete_plot.w <s> <d>; edit_plot.D <w>" blue 0 0 1
simundump xgraph /moregraphs/conc3 0 0 100 0 1 0
simundump xgraph /moregraphs/conc4 0 0 100 0 1 0
simundump xcoredraw /edit/draw 0 25 25 25 25
simundump xtree /edit/draw/tree 0 \
  /kinetics/#[],/kinetics/#[]/#[],/kinetics/#[]/#[]/#[][TYPE!=proto],/kinetics/#[]/#[]/#[][TYPE!=linkinfo]/##[] \
  "edit_elm.D <v>; drag_from_edit.w <d> <S> <x> <y> <z>" auto 0.6
simundump xtext /file/notes 0 1
addmsg /kinetics/CaM/CaM_Ca3_bind_Ca /kinetics/CaM/CaM_Ca4 REAC B A 
addmsg /kinetics/CaM/CaM_bind_Ca /kinetics/CaM/CaM REAC A B 
addmsg /kinetics/CaM/CaM_Ca3_bind_Ca /kinetics/CaM/CaM_Ca3 REAC A B 
addmsg /kinetics/CaM/CaM_Ca2_bind_Ca /kinetics/CaM/CaM_Ca3 REAC B A 
addmsg /kinetics/CaM/CaM_Ca2_bind_Ca /kinetics/CaM/CaM_Ca2 REAC A B 
addmsg /kinetics/CaM/CaM_Ca_bind_Ca /kinetics/CaM/CaM_Ca2 REAC B A 
addmsg /kinetics/CaM/CaM_bind_Ca /kinetics/CaM/CaM_Ca REAC B A 
addmsg /kinetics/CaM/CaM_Ca_bind_Ca /kinetics/CaM/CaM_Ca REAC A B 
addmsg /kinetics/CaM/CaM_Ca3 /kinetics/CaM/CaM_Ca3_bind_Ca SUBSTRATE n 
addmsg /kinetics/Ca/Ca /kinetics/CaM/CaM_Ca3_bind_Ca SUBSTRATE n 
addmsg /kinetics/CaM/CaM_Ca4 /kinetics/CaM/CaM_Ca3_bind_Ca PRODUCT n 
addmsg /kinetics/CaM/CaM /kinetics/CaM/CaM_bind_Ca SUBSTRATE n 
addmsg /kinetics/Ca/Ca /kinetics/CaM/CaM_bind_Ca SUBSTRATE n 
addmsg /kinetics/CaM/CaM_Ca /kinetics/CaM/CaM_bind_Ca PRODUCT n 
addmsg /kinetics/CaM/CaM_Ca2 /kinetics/CaM/CaM_Ca2_bind_Ca SUBSTRATE n 
addmsg /kinetics/Ca/Ca /kinetics/CaM/CaM_Ca2_bind_Ca SUBSTRATE n 
addmsg /kinetics/CaM/CaM_Ca3 /kinetics/CaM/CaM_Ca2_bind_Ca PRODUCT n 
addmsg /kinetics/CaM/CaM_Ca /kinetics/CaM/CaM_Ca_bind_Ca SUBSTRATE n 
addmsg /kinetics/Ca/Ca /kinetics/CaM/CaM_Ca_bind_Ca SUBSTRATE n 
addmsg /kinetics/CaM/CaM_Ca2 /kinetics/CaM/CaM_Ca_bind_Ca PRODUCT n 
addmsg /kinetics/CaM/CaM_Ca3_bind_Ca /kinetics/Ca/Ca REAC A B 
addmsg /kinetics/CaM/CaM_bind_Ca /kinetics/Ca/Ca REAC A B 
addmsg /kinetics/CaM/CaM_Ca2_bind_Ca /kinetics/Ca/Ca REAC A B 
addmsg /kinetics/CaM/CaM_Ca_bind_Ca /kinetics/Ca/Ca REAC A B 
addmsg /kinetics/Ca/Ca_stim /kinetics/Ca/Ca REAC B A 
addmsg /kinetics/Ca/Ca_homeostasis /kinetics/Ca/Ca REAC B A 
addmsg /kinetics/Ca/Ca_stim /kinetics/Ca/Ca_input REAC A B 
addmsg /kinetics/Ca/Ca_input /kinetics/Ca/Ca_stim SUBSTRATE n 
addmsg /kinetics/Ca/Ca /kinetics/Ca/Ca_stim PRODUCT n 
addmsg /kinetics/Ca/Ca_homeostasis /kinetics/Ca/Ca_ext REAC A B 
addmsg /kinetics/Ca/Ca_ext /kinetics/Ca/Ca_homeostasis SUBSTRATE n 
addmsg /kinetics/Ca/Ca /kinetics/Ca/Ca_homeostasis PRODUCT n 
addmsg /kinetics/Ca/Ca /graphs/conc1/Ca.Co PLOT Co *Ca *red 
addmsg /kinetics/CaM/CaM_Ca4 /graphs/conc1/CaM_Ca4.Co PLOT Co *CaM_Ca4 *blue 
enddump
// End of dump

call /kinetics/CaM/CaM/notes LOAD \
"There is a LOT of this in the cell: upto 1% of total protein mass. (Alberts et al) Say 25 uM. Meyer et al Science 256 1199-1202 1992 refer to studies saying it is comparable to CaMK levels."
call /kinetics/CaM/CaM_Ca2/notes LOAD \
"This is the intermediate where the TR2 end (the high-affinity end) has bound the Ca but the TR1 end has not."
call /kinetics/CaM/CaM_Ca/notes LOAD \
"This is the intermediate where the TR2 end (the high-affinity end) has bound the Ca but the TR1 end has not."
call /kinetics/CaM/CaM_Ca3_bind_Ca/notes LOAD \
"Use K3 = 21.5 uM here from Stemmer and Klee table 3. kb/kf =21.5 * 6e5 so kf = 7.75e-7, kb = 10"
call /kinetics/CaM/CaM_bind_Ca/notes LOAD \
"Lets use the fast rate consts here. Since the rates are so different, I am not sure whether the order is relevant. These correspond to the TR2C fragment. We use the Martin et al rates here, plus the Drabicowski binding consts. All are scaled by 3X to cell temp. kf = 2e-10 kb = 72 Stemmer & Klee: K1=.9, K2=1.1. Assume 1.0uM for both. kb/kf=3.6e11. If kb=72, kf = 2e-10 (Exactly the same !) 19 May 2006. Splitting the old CaM-TR2-bind-Ca reaction into two steps, each binding 1 Ca. This improves numerical stability and is conceptually better too. Overall rates are the same, so each kf and kb is the square root of the earlier ones. So kf = 1.125e-4, kb = 8.4853"
call /kinetics/CaM/CaM_Ca2_bind_Ca/notes LOAD \
"K3 = 21.5, K4 = 2.8. Assuming that the K4 step happens first, we get kb/kf = 2.8 uM = 1.68e6 so kf =6e-6 assuming kb = 10"
call /kinetics/CaM/CaM_Ca_bind_Ca/notes LOAD \
"Lets use the fast rate consts here. Since the rates are so different, I am not sure whether the order is relevant. These correspond to the TR2C fragment. We use the Martin et al rates here, plus the Drabicowski binding consts. All are scaled by 3X to cell temp. kf = 2e-10 kb = 72 Stemmer & Klee: K1=.9, K2=1.1. Assume 1.0uM for both. kb/kf=3.6e11. If kb=72, kf = 2e-10 (Exactly the same !) 19 May 2006. Splitting the old CaM-TR2-bind-Ca reaction into two steps, each binding 1 Ca. This improves numerical stability and is conceptually better too. Overall rates are the same, so each kf and kb is the square root of the earlier ones. So kf = 1.125e-4, kb = 8.4853"
complete_loading
