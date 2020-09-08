global_settings {
	ambient_light rgb <0.200000002980232, 0.200000002980232, 0.200000002980232>
	max_trace_level 15
}

background { color rgb <0,0,0> }

camera {
	perspective
	location <-8.02417183803928, 16.4686267408982, -8.61114427636964>
	angle 40
	up <0.456255260801853, -0.217711279176942, -0.862805271141625>
	right <0.795103757203356, 0.535106579882745, 0.285431188637993> * 1
	direction <0.399551188550027, -0.816249194232577, 0.417248248220325> }

light_source {
	<31.654660661672, 62.177061602162, -45.1474379802624>
	color rgb <1, 1, 1>
	fade_distance 113.884017101786
	fade_power 0
	parallel
	point_at <-31.654660661672, -62.177061602162, 45.1474379802624>
}

light_source {
	<-6.65815753179548, -56.2933115364977, -35.5138486853399>
	color rgb <0.300000011920929, 0.300000011920929, 0.300000011920929>
	fade_distance 113.884017101786
	fade_power 0
	parallel
	point_at <6.65815753179548, 56.2933115364977, 35.5138486853399>
}

#default {
	finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.5}
}

union {
}
union {
cylinder {
	<2.33072, -0.71345, -3.40479>, 	<2.11727635287782, -0.818984960387903, -2.96109946213837>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.11727635287782, -0.818984960387903, -2.96109946213837>, 	<1.86314, -0.94464, -2.43282>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.84821, 2.12314, -2.60775>, 	<-0.928191782019888, 1.68103037185132, -2.59294648497049>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-0.928191782019888, 1.68103037185132, -2.59294648497049>, 	<-1.02013, 1.17283, -2.57593>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<3.21108, 1.4179, -2.58285>, 	<2.97052660621131, 1.33967084580695, -2.14735139819274>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.97052660621131, 1.33967084580695, -2.14735139819274>, 	<2.68412, 1.24653, -1.62884>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-1.02013, 1.17283, -2.57593>, 	<-0.818435258163249, 1.04382237196867, -2.09381658573254>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<-0.818435258163249, 1.04382237196867, -2.09381658573254>, 	<-0.64192, 0.93092, -1.67189>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<0.82599, -0.57996, -2.46733>, 	<1.29930497500489, -0.746385787094231, -2.45158097643791>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.29930497500489, -0.746385787094231, -2.45158097643791>, 	<1.86314, -0.94464, -2.43282>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.86314, -0.94464, -2.43282>, 	<1.85937286428729, -1.54035298895886, -2.37544864312855>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.85937286428729, -1.54035298895886, -2.37544864312855>, 	<1.85621, -2.04051, -2.32728>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.86314, -0.94464, -2.43282>, 	<2.263445, -0.59879, -1.878375>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.263445, -0.59879, -1.878375>, 	<2.66375, -0.25294, -1.32393>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<4.518, -0.7196, -2.31985>, 	<4.31928067470802, -0.762547587060875, -1.85896688035013>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<4.31928067470802, -0.762547587060875, -1.85896688035013>, 	<4.08269, -0.81368, -1.31025>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.13249, 4.15116, -1.87917>, 	<0.445365, 4.255455, -1.768865>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<0.445365, 4.255455, -1.768865>, 	<1.02322, 4.35975, -1.65856>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<-2.73602, -0.06024, -1.72577>, 	<-2.84608521054057, -0.30781026101143, -1.30120797825368>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-2.84608521054057, -0.30781026101143, -1.30120797825368>, 	<-2.97713, -0.60257, -0.79572>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.65759, 1.62389, -1.72549>, 	<2.12597754073847, 1.45170728505444, -1.68139030899012>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.12597754073847, 1.45170728505444, -1.68139030899012>, 	<2.68412, 1.24653, -1.62884>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-1.91004, -5.08799, -1.69942>, 	<-1.77902960367585, -4.64401665719911, -1.51652804570377>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-1.77902960367585, -4.64401665719911, -1.51652804570377>, 	<-1.62275, -4.11441, -1.29836>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.55523, -3.19103, -1.64581>, 	<-3.09215237362105, -3.13429609396478, -1.47364297780326>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.09215237362105, -3.13429609396478, -1.47364297780326>, 	<-2.53968, -3.06661, -1.26824>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.68412, 1.24653, -1.62884>, 	<2.673935, 0.496795, -1.476385>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.673935, 0.496795, -1.476385>, 	<2.66375, -0.25294, -1.32393>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.68412, 1.24653, -1.62884>, 	<2.95988947281032, 1.56620858860883, -1.206440395445>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.95988947281032, 1.56620858860883, -1.206440395445>, 	<3.19133, 1.8345, -0.85194>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.66375, -0.25294, -1.32393>, 	<3.37322, -0.53331, -1.31709>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.37322, -0.53331, -1.31709>, 	<4.08269, -0.81368, -1.31025>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.66375, -0.25294, -1.32393>, 	<2.17255249540778, -0.377589962783792, -0.526087934799531>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.17255249540778, -0.377589962783792, -0.526087934799531>, 	<1.67304, -0.50435, 0.28526>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<4.08269, -0.81368, -1.31025>, 	<4.09721996104373, -1.39439465329711, -1.165417870539>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<4.09721996104373, -1.39439465329711, -1.165417870539>, 	<4.10942, -1.88199, -1.04381>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<4.08269, -0.81368, -1.31025>, 	<4.44134162267226, -0.51733354946914, -0.935326180100927>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<4.44134162267226, -0.51733354946914, -0.935326180100927>, 	<4.74237, -0.2686, -0.62064>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-1.62275, -4.11441, -1.29836>, 	<-2.081215, -3.59051, -1.2833>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.081215, -3.59051, -1.2833>, 	<-2.53968, -3.06661, -1.26824>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-1.62275, -4.11441, -1.29836>, 	<-0.979045, -4.011345, -1.05407>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.979045, -4.011345, -1.05407>, 	<-0.33534, -3.90828, -0.80978>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.53968, -3.06661, -1.26824>, 	<-2.33546, -2.452065, -1.016635>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.33546, -2.452065, -1.016635>, 	<-2.13124, -1.83752, -0.76503>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.80464, 1.95432, -0.93683>, 	<-3.02360444600453, 1.80279627480868, -0.511801501639966>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.02360444600453, 1.80279627480868, -0.511801501639966>, 	<-3.28448, 1.62227, -0.00542>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-4.05055, -0.8624, -0.81937>, 	<-3.56048202065384, -0.743775078651867, -0.808572634838612>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.56048202065384, -0.743775078651867, -0.808572634838612>, 	<-2.97713, -0.60257, -0.79572>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.40247, -4.71134, -0.81812>, 	<0.0660387941188283, -4.34515569212272, -0.814317074779348>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<0.0660387941188283, -4.34515569212272, -0.814317074779348>, 	<-0.33534, -3.90828, -0.80978>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.33534, -3.90828, -0.80978>, 	<-0.162685, -3.28514, -0.55104>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.162685, -3.28514, -0.55104>, 	<0.00997, -2.662, -0.2923>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.97713, -0.60257, -0.79572>, 	<-2.554185, -1.220045, -0.780375>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.554185, -1.220045, -0.780375>, 	<-2.13124, -1.83752, -0.76503>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.97713, -0.60257, -0.79572>, 	<-2.81359651569072, -0.13405754285554, -0.232189019471761>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.81359651569072, -0.13405754285554, -0.232189019471761>, 	<-2.65529, 0.31948, 0.31333>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-2.13124, -1.83752, -0.76503>, 	<-1.49724113454854, -1.74181024760962, -0.521753081657381>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-1.49724113454854, -1.74181024760962, -0.521753081657381>, 	<-0.88543, -1.64945, -0.28699>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<0.00997, -2.662, -0.2923>, 	<-0.445679244776631, -2.14673571275566, -0.289597858510427>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.445679244776631, -2.14673571275566, -0.289597858510427>, 	<-0.88543, -1.64945, -0.28699>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<0.00997, -2.662, -0.2923>, 	<0.68078, -2.509925, 0.00278>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.68078, -2.509925, 0.00278>, 	<1.35159, -2.35785, 0.29786>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.88543, -1.64945, -0.28699>, 	<-0.655563371807086, -0.719239673863844, 0.0271342871660686>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-0.655563371807086, -0.719239673863844, 0.0271342871660686>, 	<-0.40746, 0.28477, 0.36618>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<-4.35763, 1.46172, -0.23123>, 	<-3.86753099442288, 1.53504189847217, -0.128104382472748>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.86753099442288, 1.53504189847217, -0.128104382472748>, 	<-3.28448, 1.62227, -0.00542>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.14391, -2.95574, -0.17788>, 	<1.78230152696173, -2.68286781458899, 0.0392439082229729>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.78230152696173, -2.68286781458899, 0.0392439082229729>, 	<1.35159, -2.35785, 0.29786>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.28448, 1.62227, -0.00542>, 	<-2.96478850994019, 0.960322295276434, 0.15653690086709>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.96478850994019, 0.960322295276434, 0.15653690086709>, 	<-2.65529, 0.31948, 0.31333>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-3.28448, 1.62227, -0.00542>, 	<-3.21231, 2.14902, 0.532335>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.21231, 2.14902, 0.532335>, 	<-3.14014, 2.67577, 1.07009>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.67304, -0.50435, 0.28526>, 	<1.51094800748858, -1.4389821608957, 0.291613582534279>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<1.51094800748858, -1.4389821608957, 0.291613582534279>, 	<1.35159, -2.35785, 0.29786>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.67304, -0.50435, 0.28526>, 	<0.651477702013492, -0.116878122765146, 0.324993151239158>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<0.651477702013492, -0.116878122765146, 0.324993151239158>, 	<-0.40746, 0.28477, 0.36618>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<1.67304, -0.50435, 0.28526>, 	<2.24210629853049, -0.4037619288027, 1.0520738750898>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<2.24210629853049, -0.4037619288027, 1.0520738750898>, 	<2.80169, -0.30485, 1.80611>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.35159, -2.35785, 0.29786>, 	<1.34028989165074, -2.51250766376617, 0.876720143861444>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.34028989165074, -2.51250766376617, 0.876720143861444>, 	<1.3308, -2.64239, 1.36285>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-2.65529, 0.31948, 0.31333>, 	<-1.57135918467632, 0.302742418154449, 0.338814909263537>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-1.57135918467632, 0.302742418154449, 0.338814909263537>, 	<-0.40746, 0.28477, 0.36618>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<-2.65529, 0.31948, 0.31333>, 	<-2.87210966357302, 0.0477974894131194, 0.954988700814302>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-2.87210966357302, 0.0477974894131194, 0.954988700814302>, 	<-3.09606, -0.23282, 1.61775>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.40746, 0.28477, 0.36618>, 	<-0.132169698053498, 1.15543343240964, 0.579324601177586>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<-0.132169698053498, 1.15543343240964, 0.579324601177586>, 	<0.12668, 1.9741, 0.77974>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.52909, 3.63066, 0.6868>, 	<-3.35158730585558, 3.19488278541826, 0.861719675121777>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.35158730585558, 3.19488278541826, 0.861719675121777>, 	<-3.14014, 2.67577, 1.07009>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.12668, 1.9741, 0.77974>, 	<0.333557664636434, 2.52461608523938, 0.908852849432014>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.333557664636434, 2.52461608523938, 0.908852849432014>, 	<0.53054, 3.0488, 1.03179>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<4.21485, 1.15428, 0.97408>, 	<3.86945908740302, 1.11887202374495, 1.33636032104136>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.86945908740302, 1.11887202374495, 1.33636032104136>, 	<3.45799, 1.07669, 1.76795>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.14014, 2.67577, 1.07009>, 	<-2.56729747365441, 2.75490313079876, 1.21815605867196>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.56729747365441, 2.75490313079876, 1.21815605867196>, 	<-2.08658, 2.82131, 1.34241>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.14014, 2.67577, 1.07009>, 	<-3.44870031622853, 2.54377837776553, 1.56524341514023>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.44870031622853, 2.54377837776553, 1.56524341514023>, 	<-3.70774, 2.43297, 1.98093>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-5.03675, -1.11348, 1.0836>, 	<-4.83550726102491, -0.794789377146016, 1.4153748693308>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-4.83550726102491, -0.794789377146016, 1.4153748693308>, 	<-4.59579, -0.41517, 1.81058>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<4.61259, -1.38066, 1.1594>, 	<4.2739984098773, -1.38372582877766, 1.52763705946377>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<4.2739984098773, -1.38372582877766, 1.52763705946377>, 	<3.87043, -1.38738, 1.96654>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.09606, -0.23282, 1.61775>, 	<-2.81511414668222, -0.756466886145803, 1.6845846585787>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.81511414668222, -0.756466886145803, 1.6845846585787>, 	<-2.57927, -1.19605, 1.74069>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.09606, -0.23282, 1.61775>, 	<-3.845925, -0.323995, 1.714165>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.845925, -0.323995, 1.714165>, 	<-4.59579, -0.41517, 1.81058>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.09606, -0.23282, 1.61775>, 	<-2.87304820828365, 0.125400540196303, 2.04024743518303>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.87304820828365, 0.125400540196303, 2.04024743518303>, 	<-2.68589, 0.42603, 2.39482>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.71787, 1.87703, 1.62409>, 	<3.05562073364621, 1.51179811575636, 1.68973992236711>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.05562073364621, 1.51179811575636, 1.68973992236711>, 	<3.45799, 1.07669, 1.76795>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-5.14655, 0.53499, 1.74855>, 	<-4.89520327959126, 0.101371727315764, 1.77685822330408>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-4.89520327959126, 0.101371727315764, 1.77685822330408>, 	<-4.59579, -0.41517, 1.81058>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.45799, 1.07669, 1.76795>, 	<3.12984, 0.38592, 1.78703>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.12984, 0.38592, 1.78703>, 	<2.80169, -0.30485, 1.80611>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.45799, 1.07669, 1.76795>, 	<3.73610846187696, 1.17686938850748, 2.28964155736104>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.73610846187696, 1.17686938850748, 2.28964155736104>, 	<3.9697, 1.26101, 2.72781>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.80169, -0.30485, 1.80611>, 	<3.33606, -0.846115, 1.886325>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.33606, -0.846115, 1.886325>, 	<3.87043, -1.38738, 1.96654>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.80169, -0.30485, 1.80611>, 	<2.331015, -0.33359, 2.41029>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.331015, -0.33359, 2.41029>, 	<1.86034, -0.36233, 3.01447>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-4.59579, -0.41517, 1.81058>, 	<-4.69352295455317, -0.643046848868597, 2.35490656616164>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-4.69352295455317, -0.643046848868597, 2.35490656616164>, 	<-4.77557, -0.83435, 2.81187>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.87043, -1.38738, 1.96654>, 	<3.63472450416095, -1.93635380481105, 1.9984066393738>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.63472450416095, -1.93635380481105, 1.9984066393738>, 	<3.43684, -2.39724, 2.02516>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.87043, -1.38738, 1.96654>, 	<4.16494072766129, -1.29448531192909, 2.48067285791272>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<4.16494072766129, -1.29448531192909, 2.48067285791272>, 	<4.41231, -1.21646, 2.91251>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.06471, 0.39181, 2.93114>, 	<1.42778790764089, 0.0476656448747533, 2.96916682408119>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.42778790764089, 0.0476656448747533, 2.96916682408119>, 	<1.86034, -0.36233, 3.01447>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.86034, -0.36233, 3.01447>, 	<1.59478826909592, -0.89655894273643, 3.06748902778152>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.59478826909592, -0.89655894273643, 3.06748902778152>, 	<1.3718, -1.34516, 3.11201>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.86034, -0.36233, 3.01447>, 	<2.17588501769541, -0.267070819092158, 3.51550057289726>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.17588501769541, -0.267070819092158, 3.51550057289726>, 	<2.44092, -0.18706, 3.93633>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
sphere {
	<-3.52909, 3.63066, 0.6868>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-3.14014, 2.67577, 1.07009>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-3.70774, 2.43297, 1.98093>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.08658, 2.82131, 1.34241>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.65759, 1.62389, -1.72549>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-0.13249, 4.15116, -1.87917>, 0.2432
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581,0> }
}
sphere {
	<4.21485, 1.15428, 0.97408>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.02322, 4.35975, -1.65856>, 0.2432
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581,0> }
}
sphere {
	<-3.28448, 1.62227, -0.00542>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<3.21108, 1.4179, -2.58285>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-4.35763, 1.46172, -0.23123>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.68412, 1.24653, -1.62884>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<2.71787, 1.87703, 1.62409>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<3.19133, 1.8345, -0.85194>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.80464, 1.95432, -0.93683>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<3.45799, 1.07669, 1.76795>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<3.9697, 1.26101, 2.72781>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.73602, -0.06024, -1.72577>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-0.40746, 0.28477, 0.36618>, 0.328
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186,0> }
}
sphere {
	<0.82599, -0.57996, -2.46733>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.65529, 0.31948, 0.31333>, 0.248
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1,0> }
}
sphere {
	<-2.97713, -0.60257, -0.79572>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<2.66375, -0.25294, -1.32393>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<2.33072, -0.71345, -3.40479>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<4.61259, -1.38066, 1.1594>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.86314, -0.94464, -2.43282>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<2.80169, -0.30485, 1.80611>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<1.67304, -0.50435, 0.28526>, 0.288
	pigment { rgbt <1, 0.5, 0,0> }
}
sphere {
	<-4.05055, -0.8624, -0.81937>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<3.87043, -1.38738, 1.96654>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-2.68589, 0.42603, 2.39482>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-5.14655, 0.53499, 1.74855>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<4.41231, -1.21646, 2.91251>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-3.09606, -0.23282, 1.61775>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<4.518, -0.7196, -2.31985>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<0.12668, 1.9741, 0.77974>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<4.08269, -0.81368, -1.31025>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<1.06471, 0.39181, 2.93114>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.86034, -0.36233, 3.01447>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<4.74237, -0.2686, -0.62064>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.44092, -0.18706, 3.93633>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.13124, -1.83752, -0.76503>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-0.88543, -1.64945, -0.28699>, 0.248
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1,0> }
}
sphere {
	<-4.59579, -0.41517, 1.81058>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<1.85621, -2.04051, -2.32728>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.57927, -1.19605, 1.74069>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<3.43684, -2.39724, 2.02516>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-5.03675, -1.11348, 1.0836>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.35159, -2.35785, 0.29786>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-4.77557, -0.83435, 2.81187>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<4.10942, -1.88199, -1.04381>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<0.00997, -2.662, -0.2923>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-2.53968, -3.06661, -1.26824>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<1.3718, -1.34516, 3.11201>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-3.55523, -3.19103, -1.64581>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.3308, -2.64239, 1.36285>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.14391, -2.95574, -0.17788>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-0.33534, -3.90828, -0.80978>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-1.62275, -4.11441, -1.29836>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<0.40247, -4.71134, -0.81812>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-1.91004, -5.08799, -1.69942>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<0.53054, 3.0488, 1.03179>, 0.2432
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581,0> }
}
sphere {
	<-1.02013, 1.17283, -2.57593>, 0.2432
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581,0> }
}
sphere {
	<-0.84821, 2.12314, -2.60775>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-0.64192, 0.93092, -1.67189>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
}
merge {
}
