global_settings {
	ambient_light rgb <0.200000002980232, 0.200000002980232, 0.200000002980232>
	max_trace_level 15
}

background { color rgb <0,0,0> }

camera {
	perspective
	location <5.12768847230217, -13.0250717369394, -8.44330909611656>
	angle 40
	up <-0.0547223973724771, 0.517789346777145, -0.853756201494259>
	right <-0.949400136273418, -0.291811697747068, -0.116126286007906> * 1
	direction <-0.309265000395817, 0.804201545274611, 0.507557912073198> }

light_source {
	<-24.6028766812267, -34.0015456183707, -60.3331626094533>
	color rgb <1, 1, 1>
	fade_distance 100.715721808572
	fade_power 0
	parallel
	point_at <24.6028766812267, 34.0015456183707, 60.3331626094533>
}

light_source {
	<28.531854657851, 50.2572796825395, -12.6372721914616>
	color rgb <0.300000011920929, 0.300000011920929, 0.300000011920929>
	fade_distance 100.715721808572
	fade_power 0
	parallel
	point_at <-28.531854657851, -50.2572796825395, 12.6372721914616>
}

#default {
	finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.5}
}

union {
}
union {
cylinder {
	<-0.69487285, -2.18311611, -3.27277416>, 	<-0.53395241077008, -2.26664965175599, -2.85789615111754>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-0.53395241077008, -2.26664965175599, -2.85789615111754>, 	<-0.34914799, -2.36258133, -2.38144151>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<-2.12184913, 1.04499984, -2.41770806>, 	<-2.50196363745644, 0.911194205040454, -2.13172828909695>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-2.50196363745644, 0.911194205040454, -2.13172828909695>, 	<-2.95592546, 0.75139328, -1.79018937>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.34914799, -2.36258133, -2.38144151>, 	<-0.11614571, -1.592285885, -2.20865644>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<-0.11614571, -1.592285885, -2.20865644>, 	<0.11685657, -0.82199044, -2.03587137>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<-3.424951, -0.11461877, -2.24627638>, 	<-3.21118112119389, 0.280087435067484, -2.03840355317539>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.21118112119389, 0.280087435067484, -2.03840355317539>, 	<-2.95592546, 0.75139328, -1.79018937>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.11685657, -0.82199044, -2.03587137>, 	<0.00460298394979224, -0.850523133276259, -1.03057299819318>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<0.00460298394979224, -0.850523133276259, -1.03057299819318>, 	<-0.11705731, -0.88144683, 0.0589681>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<-2.95592546, 0.75139328, -1.79018937>, 	<-3.34958797519334, 1.19127265166371, -1.75628161675173>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.34958797519334, 1.19127265166371, -1.75628161675173>, 	<-3.67933575, 1.55973355, -1.7278791>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-2.95592546, 0.75139328, -1.79018937>, 	<-2.63471938676144, 0.514287489998184, -0.95542140352388>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.63471938676144, 0.514287489998184, -0.95542140352388>, 	<-2.30795928, 0.27308186, -0.10621931>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<2.47961778, -0.22805005, -1.63105341>, 	<2.54760573632717, -0.151619199146972, -1.14776565637381>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.54760573632717, -0.151619199146972, -1.14776565637381>, 	<2.62880606, -0.06033524, -0.57055869>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.00477094, -2.45857026, -1.57612531>, 	<2.21202836685001, -2.47184195894365, -1.1300384429631>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.21202836685001, -2.47184195894365, -1.1300384429631>, 	<2.45972115, -2.48770293, -0.59692125>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.54399273, -2.58737731, -0.68353774>, 	<3.04950487880292, -2.54192027608457, -0.644035826634039>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.04950487880292, -2.54192027608457, -0.644035826634039>, 	<2.45972115, -2.48770293, -0.59692125>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.45972115, -2.48770293, -0.59692125>, 	<2.24874738361894, -2.95268183921413, -0.301587661648098>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.24874738361894, -2.95268183921413, -0.301587661648098>, 	<2.07210958, -3.34198547, -0.05431959>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.45972115, -2.48770293, -0.59692125>, 	<2.27887844366883, -1.86001567983372, -0.221197771523344>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.27887844366883, -1.86001567983372, -0.221197771523344>, 	<2.10379531, -1.25231934, 0.14255947>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<2.62880606, -0.06033524, -0.57055869>, 	<3.20788662305422, 0.00604770142862372, -0.465700516028258>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.20788662305422, 0.00604770142862372, -0.465700516028258>, 	<3.69309868, 0.06167002, -0.37783977>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.62880606, -0.06033524, -0.57055869>, 	<2.25806755, 0.57255079, -0.386627335>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.25806755, 0.57255079, -0.386627335>, 	<1.88732904, 1.20543682, -0.20269598>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.62880606, -0.06033524, -0.57055869>, 	<2.36205795733623, -0.665959976408496, -0.208236744868255>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.36205795733623, -0.665959976408496, -0.208236744868255>, 	<2.10379531, -1.25231934, 0.14255947>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<3.52977925, 2.55859074, -0.42706967>, 	<3.05117095784357, 2.51628750474437, -0.333856023895179>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.05117095784357, 2.51628750474437, -0.333856023895179>, 	<2.4786858, 2.46568668, -0.22235894>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.4786858, 2.46568668, -0.22235894>, 	<2.18300742, 1.83556175, -0.21252746>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.18300742, 1.83556175, -0.21252746>, 	<1.88732904, 1.20543682, -0.20269598>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.4786858, 2.46568668, -0.22235894>, 	<2.090372875, 3.020520975, -0.106140755>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.090372875, 3.020520975, -0.106140755>, 	<1.70205995, 3.57535527, 0.01007743>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.88732904, 1.20543682, -0.20269598>, 	<1.22163247511895, 1.13826225749144, -0.0671319443185073>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.22163247511895, 1.13826225749144, -0.0671319443185073>, 	<0.5793391, 1.07344928, 0.06366621>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-2.30795928, 0.27308186, -0.10621931>, 	<-1.2301626923865, -0.294879235212253, -0.0249566933915938>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<-1.2301626923865, -0.294879235212253, -0.0249566933915938>, 	<-0.11705731, -0.88144683, 0.0589681>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<-2.30795928, 0.27308186, -0.10621931>, 	<-1.98335465600131, 1.09322065835942, 0.264233610182045>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<-1.98335465600131, 1.09322065835942, 0.264233610182045>, 	<-1.6641789, 1.899643, 0.62849087>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.30795928, 0.27308186, -0.10621931>, 	<-3.08657151173996, 0.121603079272596, 0.401072265604342>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<-3.08657151173996, 0.121603079272596, 0.401072265604342>, 	<-3.85195267, -0.0273016, 0.89974336>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.13906704, 4.55889235, -0.0058556>, 	<1.94005151807819, 4.11098395266565, 0.00140039276031382>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.94005151807819, 4.11098395266565, 0.00140039276031382>, 	<1.70205995, 3.57535527, 0.01007743>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.70205995, 3.57535527, 0.01007743>, 	<1.021424205, 3.498268825, 0.138794785>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.021424205, 3.498268825, 0.138794785>, 	<0.34078846, 3.42118238, 0.26751214>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.11705731, -0.88144683, 0.0589681>, 	<0.244563870624681, 0.133681664554986, 0.0614077106305382>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<0.244563870624681, 0.133681664554986, 0.0614077106305382>, 	<0.5793391, 1.07344928, 0.06366621>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-0.11705731, -0.88144683, 0.0589681>, 	<1.03279548972092, -1.07346713342879, 0.10224777104229>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<1.03279548972092, -1.07346713342879, 0.10224777104229>, 	<2.10379531, -1.25231934, 0.14255947>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-0.11705731, -0.88144683, 0.0589681>, 	<-0.269837648798794, -0.965440481738265, 1.15395941600081>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<-0.269837648798794, -0.965440481738265, 1.15395941600081>, 	<-0.4149017, -1.04519197, 2.19364736>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.5793391, 1.07344928, 0.06366621>, 	<0.191902775619962, 1.60018453445452, 0.179895353503272>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<0.191902775619962, 1.60018453445452, 0.179895353503272>, 	<-0.20953457, 2.14595474, 0.30032474>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.10379531, -1.25231934, 0.14255947>, 	<2.36184955744489, -1.28054654359285, 0.823167445247086>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<2.36184955744489, -1.28054654359285, 0.823167445247086>, 	<2.62840602, -1.30970376, 1.52619968>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.34078846, 3.42118238, 0.26751214>, 	<0.065626945, 2.78356856, 0.28391844>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.065626945, 2.78356856, 0.28391844>, 	<-0.20953457, 2.14595474, 0.30032474>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.34078846, 3.42118238, 0.26751214>, 	<0.000778862824584805, 3.88695245365828, 0.367148324389787>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.000778862824584805, 3.88695245365828, 0.367148324389787>, 	<-0.28345471, 4.27631644, 0.45043994>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-0.20953457, 2.14595474, 0.30032474>, 	<-0.936856735, 2.02279887, 0.464407805>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.936856735, 2.02279887, 0.464407805>, 	<-1.6641789, 1.899643, 0.62849087>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.26840551, 2.74654963, 0.32032961>, 	<-1.99302207062388, 2.36056189503441, 0.460777756865425>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-1.99302207062388, 2.36056189503441, 0.460777756865425>, 	<-1.6641789, 1.899643, 0.62849087>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-4.38840515, -0.87507903, 0.48484844>, 	<-4.14389706319311, -0.488673066423584, 0.673952166606121>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-4.14389706319311, -0.488673066423584, 0.673952166606121>, 	<-3.85195267, -0.0273016, 0.89974336>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-1.6641789, 1.899643, 0.62849087>, 	<-1.71514316412599, 1.83996454209091, 1.21563609364491>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-1.71514316412599, 1.83996454209091, 1.21563609364491>, 	<-1.75784807, 1.78995768, 1.70762753>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-4.50670963, 0.83997556, 0.89817174>, 	<-4.20825250887981, 0.444645597586324, 0.898888129758873>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-4.20825250887981, 0.444645597586324, 0.898888129758873>, 	<-3.85195267, -0.0273016, 0.89974336>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.85195267, -0.0273016, 0.89974336>, 	<-3.69975318625712, -0.158828148188227, 1.45483550484822>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.69975318625712, -0.158828148188227, 1.45483550484822>, 	<-3.57230091, -0.26896885, 1.91967123>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.71628435, -1.40166085, 1.51688219>, 	<3.22017287377667, -1.35972512390727, 1.52113129910634>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.22017287377667, -1.35972512390727, 1.52113129910634>, 	<2.62840602, -1.30970376, 1.52619968>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.62840602, -1.30970376, 1.52619968>, 	<2.39554102421963, -1.77278659510196, 1.80518866328453>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.39554102421963, -1.77278659510196, 1.80518866328453>, 	<2.20065294, -2.16034732, 2.03867858>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.62840602, -1.30970376, 1.52619968>, 	<2.47481081291144, -0.822808739058436, 1.82015536622942>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.47481081291144, -0.822808739058436, 1.82015536622942>, 	<2.34624449, -0.41525499, 2.06620993>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-0.4149017, -1.04519197, 2.19364736>, 	<-0.576509455962133, -1.4039775649776, 2.67323432634079>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.576509455962133, -1.4039775649776, 2.67323432634079>, 	<-0.73061464, -1.74610669, 3.13055667>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
sphere {
	<-0.34914799, -2.36258133, -2.38144151>, 0.2432
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581,0> }
}
sphere {
	<0.11685657, -0.82199044, -2.03587137>, 0.2432
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581,0> }
}
sphere {
	<2.45972115, -2.48770293, -0.59692125>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<2.47961778, -0.22805005, -1.63105341>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.10379531, -1.25231934, 0.14255947>, 0.248
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1,0> }
}
sphere {
	<2.62880606, -0.06033524, -0.57055869>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-2.95592546, 0.75139328, -1.79018937>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-3.85195267, -0.0273016, 0.89974336>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-2.30795928, 0.27308186, -0.10621931>, 0.288
	pigment { rgbt <1, 0.5, 0,0> }
}
sphere {
	<3.69309868, 0.06167002, -0.37783977>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.62840602, -1.30970376, 1.52619968>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-0.4149017, -1.04519197, 2.19364736>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<1.88732904, 1.20543682, -0.20269598>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<0.5793391, 1.07344928, 0.06366621>, 0.248
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1,0> }
}
sphere {
	<-1.6641789, 1.899643, 0.62849087>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-0.20953457, 2.14595474, 0.30032474>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<2.4786858, 2.46568668, -0.22235894>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<3.52977925, 2.55859074, -0.42706967>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-1.75784807, 1.78995768, 1.70762753>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.26840551, 2.74654963, 0.32032961>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<0.34078846, 3.42118238, 0.26751214>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<1.70205995, 3.57535527, 0.01007743>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-0.28345471, 4.27631644, 0.45043994>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.13906704, 4.55889235, -0.0058556>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-0.73061464, -1.74610669, 3.13055667>, 0.2432
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581,0> }
}
sphere {
	<-3.67933575, 1.55973355, -1.7278791>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-3.424951, -0.11461877, -2.24627638>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.12184913, 1.04499984, -2.41770806>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-3.57230091, -0.26896885, 1.91967123>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-4.38840515, -0.87507903, 0.48484844>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-4.50670963, 0.83997556, 0.89817174>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.07210958, -3.34198547, -0.05431959>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<3.54399273, -2.58737731, -0.68353774>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.00477094, -2.45857026, -1.57612531>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.34624449, -0.41525499, 2.06620993>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<3.71628435, -1.40166085, 1.51688219>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.20065294, -2.16034732, 2.03867858>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-0.11705731, -0.88144683, 0.0589681>, 0.328
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186,0> }
}
sphere {
	<-0.69487285, -2.18311611, -3.27277416>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
}
merge {
}
