global_settings {
	ambient_light rgb <0.200000002980232, 0.200000002980232, 0.200000002980232>
	max_trace_level 15
}

background { color rgb <0,0,0> }

camera {
	perspective
	location <-2.08831718888988, -15.0039517664006, 9.20445108618485>
	angle 40
	up <0.156507220406495, 0.491718688975785, 0.856573535005938>
	right <0.981762269392441, -0.172220250155225, -0.0805172766173495> * 1
	direction <0.107927458776885, 0.853553112787018, -0.509704568641048> }

light_source {
	<43.8413284660934, -36.0476811280554, 58.2049358848482>
	color rgb <1, 1, 1>
	fade_distance 111.408498005644
	fade_power 0
	parallel
	point_at <-43.8413284660934, 36.0476811280554, -58.2049358848482>
}

light_source {
	<-34.6419736461559, 50.6215613151645, 22.7921177587512>
	color rgb <0.300000011920929, 0.300000011920929, 0.300000011920929>
	fade_distance 111.408498005644
	fade_power 0
	parallel
	point_at <34.6419736461559, -50.6215613151645, -22.7921177587512>
}

#default {
	finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.5}
}

union {
}
union {
cylinder {
	<2.904345, -1.968492, -3.00088>, 	<2.58080285082086, -1.78759187101034, -2.66016321878571>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.58080285082086, -1.78759187101034, -2.66016321878571>, 	<2.19557, -1.572199, -2.254481>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.45853, -0.921121, -2.745631>, 	<1.79486642143612, -1.21823003695293, -2.52150254918546>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.79486642143612, -1.21823003695293, -2.52150254918546>, 	<2.19557, -1.572199, -2.254481>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.453452, 0.788738, -2.563694>, 	<3.68790728706034, 0.485088105944692, -2.24011423371594>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.68790728706034, 0.485088105944692, -2.24011423371594>, 	<3.96721, 0.123355, -1.854639>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<4.692695, -0.469255, -2.436925>, 	<4.36152440240453, -0.198739403066842, -2.17112211973166>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<4.36152440240453, -0.198739403066842, -2.17112211973166>, 	<3.96721, 0.123355, -1.854639>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.673779, -1.236077, -2.419907>, 	<-0.649800486022483, -0.830584662975413, -2.23041664445025>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-0.649800486022483, -0.830584662975413, -2.23041664445025>, 	<-0.622227, -0.364299, -2.012517>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<2.19557, -1.572199, -2.254481>, 	<1.89369750559935, -2.0315057002967, -2.02071324559864>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.89369750559935, -2.0315057002967, -2.02071324559864>, 	<1.640345, -2.416988, -1.824519>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.19557, -1.572199, -2.254481>, 	<2.591985, -1.2035145, -1.712474>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.591985, -1.2035145, -1.712474>, 	<2.9884, -0.83483, -1.170467>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-5.157072, -0.323684, -2.080738>, 	<-4.99585367282305, -0.273929685640309, -1.60769435743921>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-4.99585367282305, -0.273929685640309, -1.60769435743921>, 	<-4.803819, -0.214665, -1.04423>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.205324, 1.245143, -2.063176>, 	<1.28514036140271, 1.42414886984038, -1.5986359198187>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.28514036140271, 1.42414886984038, -1.5986359198187>, 	<1.380154, 1.637238, -1.045646>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.622227, -0.364299, -2.012517>, 	<-0.517582173061353, -0.493944520520098, -0.987030747146845>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<-0.517582173061353, -0.493944520520098, -0.987030747146845>, 	<-0.404396, -0.634172, 0.122158>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<3.96721, 0.123355, -1.854639>, 	<3.477805, -0.3557375, -1.512553>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.477805, -0.3557375, -1.512553>, 	<2.9884, -0.83483, -1.170467>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.96721, 0.123355, -1.854639>, 	<4.27707869174395, 0.459658026377219, -1.47149258178238>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<4.27707869174395, 0.459658026377219, -1.47149258178238>, 	<4.537085, 0.741845, -1.15>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-2.723029, 0.148521, -1.636366>, 	<-2.99397499005115, -0.160608725590299, -1.34874902222606>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-2.99397499005115, -0.160608725590299, -1.34874902222606>, 	<-3.316767, -0.528891, -1.006096>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.128001, -1.520549, -1.435864>, 	<-3.21412547472708, -1.06810510460092, -1.23978237908042>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.21412547472708, -1.06810510460092, -1.23978237908042>, 	<-3.316767, -0.528891, -1.006096>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.48999, 4.273662, -1.246163>, 	<0.158924999311699, 3.96563999884311, -1.03941392162787>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<0.158924999311699, 3.96563999884311, -1.03941392162787>, 	<-0.23606, 3.598147, -0.792747>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.9884, -0.83483, -1.170467>, 	<3.3695845, -1.356348, -0.7607665>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.3695845, -1.356348, -0.7607665>, 	<3.750769, -1.877866, -0.351066>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.9884, -0.83483, -1.170467>, 	<2.36371380236676, -0.383154573337373, -0.627859752111301>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.36371380236676, -0.383154573337373, -0.627859752111301>, 	<1.728425, 0.076187, -0.076043>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<2.200108, 2.369213, -1.066183>, 	<1.82593250948445, 2.03518559764558, -1.05681120347395>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.82593250948445, 2.03518559764558, -1.05681120347395>, 	<1.380154, 1.637238, -1.045646>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.380154, 1.637238, -1.045646>, 	<0.733626, 1.945156, -0.828516>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.733626, 1.945156, -0.828516>, 	<0.087098, 2.253074, -0.611386>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.380154, 1.637238, -1.045646>, 	<1.55279986936441, 0.86338945384719, -0.564991702417536>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.55279986936441, 0.86338945384719, -0.564991702417536>, 	<1.728425, 0.076187, -0.076043>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<-4.803819, -0.214665, -1.04423>, 	<-4.060293, -0.371778, -1.025163>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-4.060293, -0.371778, -1.025163>, 	<-3.316767, -0.528891, -1.006096>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-4.803819, -0.214665, -1.04423>, 	<-4.92812283007237, 0.346599694573325, -0.878152168109508>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-4.92812283007237, 0.346599694573325, -0.878152168109508>, 	<-5.032484, 0.817818, -0.738719>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-4.803819, -0.214665, -1.04423>, 	<-5.13083081304322, -0.584568402378816, -0.705532884242823>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-5.13083081304322, -0.584568402378816, -0.705532884242823>, 	<-5.40541, -0.895162, -0.421142>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<4.336951, -2.506362, -1.042229>, 	<4.06937228058514, -2.21946790424584, -0.726728847407568>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<4.06937228058514, -2.21946790424584, -0.726728847407568>, 	<3.750769, -1.877866, -0.351066>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.316767, -0.528891, -1.006096>, 	<-2.99048737731928, -0.529307129566754, -0.327207794907634>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.99048737731928, -0.529307129566754, -0.327207794907634>, 	<-2.674604, -0.52971, 0.330049>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-0.23606, 3.598147, -0.792747>, 	<-0.074481, 2.9256105, -0.7020665>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.074481, 2.9256105, -0.7020665>, 	<0.087098, 2.253074, -0.611386>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.23606, 3.598147, -0.792747>, 	<-0.864433, 3.825519, -0.5996695>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.864433, 3.825519, -0.5996695>, 	<-1.492806, 4.052891, -0.406592>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.087098, 2.253074, -0.611386>, 	<-0.364419788934625, 1.83011524801749, -0.321284407035145>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.364419788934625, 1.83011524801749, -0.321284407035145>, 	<-0.800079, 1.422012, -0.041372>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-1.764558, 5.100066, -0.552391>, 	<-1.64063106161497, 4.62252327095533, -0.48590234972475>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-1.64063106161497, 4.62252327095533, -0.48590234972475>, 	<-1.492806, 4.052891, -0.406592>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-1.492806, 4.052891, -0.406592>, 	<-1.9494635, 3.608388, -0.124561>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-1.9494635, 3.608388, -0.124561>, 	<-2.406121, 3.163885, 0.15747>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.750769, -1.877866, -0.351066>, 	<3.37809224000593, -2.23879113021471, -0.0534774627770804>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.37809224000593, -2.23879113021471, -0.0534774627770804>, 	<3.065234, -2.541784, 0.196345>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.750769, -1.877866, -0.351066>, 	<4.13391559101126, -1.63545028197667, 0.037859448818008>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<4.13391559101126, -1.63545028197667, 0.037859448818008>, 	<4.455487, -1.431993, 0.364281>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-4.05502, -3.048759, -0.13346>, 	<-3.7078535659493, -3.00377598959379, 0.226205568391941>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.7078535659493, -3.00377598959379, 0.226205568391941>, 	<-3.294284, -2.950189, 0.654665>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.728425, 0.076187, -0.076043>, 	<0.680916392316358, -0.272697021230823, 0.0213009653639502>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<0.680916392316358, -0.272697021230823, 0.0213009653639502>, 	<-0.404396, -0.634172, 0.122158>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<1.728425, 0.076187, -0.076043>, 	<2.15172828314381, 0.375326153698926, 0.732076865989208>, 0.05
	pigment { rgbt <1, 0.5, 0, 0> }
}
cylinder {
	<2.15172828314381, 0.375326153698926, 0.732076865989208>, 	<2.567977, 0.66948, 1.526729>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.800079, 1.422012, -0.041372>, 	<-0.609773294062011, 0.433080083141308, 0.0372785664687119>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-0.609773294062011, 0.433080083141308, 0.0372785664687119>, 	<-0.404396, -0.634172, 0.122158>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<-0.800079, 1.422012, -0.041372>, 	<-1.39679691459309, 1.62826719407269, 0.142526171645841>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-1.39679691459309, 1.62826719407269, 0.142526171645841>, 	<-2.015293, 1.84205, 0.333136>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.404396, -0.634172, 0.122158>, 	<-1.57929157977131, -0.580110018871367, 0.229747356118135>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<-1.57929157977131, -0.580110018871367, 0.229747356118135>, 	<-2.674604, -0.52971, 0.330049>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-0.404396, -0.634172, 0.122158>, 	<-0.180703794261403, -1.54066294181291, 0.232123080543017>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<-0.180703794261403, -1.54066294181291, 0.232123080543017>, 	<0.029664, -2.393158, 0.335538>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-0.404396, -0.634172, 0.122158>, 	<-0.39345046480913, -0.592485789515112, 1.0088693831536>, 0.05
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186, 0> }
}
cylinder {
	<-0.39345046480913, -0.592485789515112, 1.0088693831536>, 	<-0.384379, -0.557937, 1.74376>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-2.406121, 3.163885, 0.15747>, 	<-2.210707, 2.5029675, 0.245303>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.210707, 2.5029675, 0.245303>, 	<-2.015293, 1.84205, 0.333136>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.406121, 3.163885, 0.15747>, 	<-2.94785524080685, 3.34050261709006, 0.322483133564119>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.94785524080685, 3.34050261709006, 0.322483133564119>, 	<-3.40192, 3.488538, 0.460792>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-2.321065, -3.247014, 0.240905>, 	<-2.76514446135289, -3.1115728648433, 0.429703531398763>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-2.76514446135289, -3.1115728648433, 0.429703531398763>, 	<-3.294284, -2.950189, 0.654665>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.674604, -0.52971, 0.330049>, 	<-2.75406065110968, 0.11416589850219, 0.657122871588636>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-2.75406065110968, 0.11416589850219, 0.657122871588636>, 	<-2.836142, 0.779311, 0.995001>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.674604, -0.52971, 0.330049>, 	<-2.95214793703003, -1.03395886010151, 0.77652900348093>, 0.05
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1, 0> }
}
cylinder {
	<-2.95214793703003, -1.03395886010151, 0.77652900348093>, 	<-3.238836, -1.554821, 1.237719>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.015293, 1.84205, 0.333136>, 	<-2.4257175, 1.3106805, 0.6640685>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.4257175, 1.3106805, 0.6640685>, 	<-2.836142, 0.779311, 0.995001>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.029664, -2.393158, 0.335538>, 	<0.199911239753415, -2.95967023426456, 0.423055005213556>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<0.199911239753415, -2.95967023426456, 0.423055005213556>, 	<0.361959, -3.498898, 0.506357>, 0.05
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581, 0> }
}
cylinder {
	<3.95127, 2.103295, 0.60632>, 	<3.95471409001126, 1.72217672778952, 0.933982307421384>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.95471409001126, 1.72217672778952, 0.933982307421384>, 	<3.958816, 1.268265, 1.324228>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.294284, -2.950189, 0.654665>, 	<-3.26656, -2.252505, 0.946192>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.26656, -2.252505, 0.946192>, 	<-3.238836, -1.554821, 1.237719>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.294284, -2.950189, 0.654665>, 	<-3.43488721075254, -3.33673053926731, 1.08867826008149>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.43488721075254, -3.33673053926731, 1.08867826008149>, 	<-3.552917, -3.661214, 1.453012>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<4.696633, 0.528244, 0.986588>, 	<4.35997342194601, 0.865909244377793, 1.14065024027658>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<4.35997342194601, 0.865909244377793, 1.14065024027658>, 	<3.958816, 1.268265, 1.324228>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.836142, 0.779311, 0.995001>, 	<-3.41213727127358, 0.942302051878826, 1.03129871283475>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.41213727127358, 0.942302051878826, 1.03129871283475>, 	<-3.895929, 1.079202, 1.061786>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-2.836142, 0.779311, 0.995001>, 	<-2.62749299909277, 0.712687991162995, 1.55260846078599>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.62749299909277, 0.712687991162995, 1.55260846078599>, 	<-2.452279, 0.656741, 2.020862>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.238836, -1.554821, 1.237719>, 	<-3.79054964618267, -1.38700330226197, 1.41075307097392>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-3.79054964618267, -1.38700330226197, 1.41075307097392>, 	<-4.254293, -1.245944, 1.556197>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<-3.238836, -1.554821, 1.237719>, 	<-2.88899048893167, -1.55060329519898, 1.72233915315238>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<-2.88899048893167, -1.55060329519898, 1.72233915315238>, 	<-2.595334, -1.547063, 2.129124>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.958816, 1.268265, 1.324228>, 	<3.2633965, 0.9688725, 1.4254785>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.2633965, 0.9688725, 1.4254785>, 	<2.567977, 0.66948, 1.526729>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<3.958816, 1.268265, 1.324228>, 	<4.15404771326287, 1.48710905207589, 1.84717847471192>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<4.15404771326287, 1.48710905207589, 1.84717847471192>, 	<4.31802, 1.670913, 2.286397>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.63082, 2.66104, 1.516873>, 	<1.64551019490583, 2.24454093603006, 1.79725956743331>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.64551019490583, 2.24454093603006, 1.79725956743331>, 	<1.663008, 1.74844, 2.131234>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.567977, 0.66948, 1.526729>, 	<2.1154925, 1.20896, 1.8289815>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.1154925, 1.20896, 1.8289815>, 	<1.663008, 1.74844, 2.131234>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.567977, 0.66948, 1.526729>, 	<2.595344, 0.084074, 2.0193565>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.595344, 0.084074, 2.0193565>, 	<2.622711, -0.501332, 2.511984>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.663008, 1.74844, 2.131234>, 	<1.10613411241227, 1.54430055138649, 2.20427752328208>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.10613411241227, 1.54430055138649, 2.20427752328208>, 	<0.638719, 1.372955, 2.265587>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<1.663008, 1.74844, 2.131234>, 	<1.88051955176094, 1.90172079523954, 2.6682469983874>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<1.88051955176094, 1.90172079523954, 2.6682469983874>, 	<2.06319, 2.030449, 3.119241>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<3.310396, -1.294196, 2.191361>, 	<2.99663196550117, -0.932442860919059, 2.33764885858804>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.99663196550117, -0.932442860919059, 2.33764885858804>, 	<2.622711, -0.501332, 2.511984>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.622711, -0.501332, 2.511984>, 	<2.07866885698276, -0.739486824035008, 2.58128599260284>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.07866885698276, -0.739486824035008, 2.58128599260284>, 	<1.621977, -0.939404, 2.639461>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
cylinder {
	<2.622711, -0.501332, 2.511984>, 	<2.81222149094961, -0.30074842292376, 3.04419816357208>, 0.05
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464, 0> }
}
cylinder {
	<2.81222149094961, -0.30074842292376, 3.04419816357208>, 	<2.971385, -0.132285, 3.491187>, 0.05
	pigment { rgbt <0.75, 0.75, 0.75, 0> }
}
sphere {
	<-0.800079, 1.422012, -0.041372>, 0.248
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1,0> }
}
sphere {
	<-2.015293, 1.84205, 0.333136>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<0.087098, 2.253074, -0.611386>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-2.406121, 3.163885, 0.15747>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-1.492806, 4.052891, -0.406592>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-0.23606, 3.598147, -0.792747>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-2.836142, 0.779311, 0.995001>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-3.40192, 3.488538, 0.460792>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-1.764558, 5.100066, -0.552391>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<0.48999, 4.273662, -1.246163>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.380154, 1.637238, -1.045646>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-2.674604, -0.52971, 0.330049>, 0.248
	pigment { rgbt <0.0500000007450581, 0.0500000007450581, 1,0> }
}
sphere {
	<-3.895929, 1.079202, 1.061786>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.452279, 0.656741, 2.020862>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.205324, 1.245143, -2.063176>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.728425, 0.076187, -0.076043>, 0.288
	pigment { rgbt <1, 0.5, 0,0> }
}
sphere {
	<2.200108, 2.369213, -1.066183>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-0.404396, -0.634172, 0.122158>, 0.328
	pigment { rgbt <0.140000000596046, 0.560000002384186, 0.560000002384186,0> }
}
sphere {
	<-3.316767, -0.528891, -1.006096>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<2.9884, -0.83483, -1.170467>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-3.238836, -1.554821, 1.237719>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<2.567977, 0.66948, 1.526729>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<2.19557, -1.572199, -2.254481>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<3.750769, -1.877866, -0.351066>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<3.96721, 0.123355, -1.854639>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<1.663008, 1.74844, 2.131234>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<3.958816, 1.268265, 1.324228>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<2.622711, -0.501332, 2.511984>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-0.384379, -0.557937, 1.74376>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.723029, 0.148521, -1.636366>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-4.803819, -0.214665, -1.04423>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-3.128001, -1.520549, -1.435864>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.595334, -1.547063, 2.129124>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-3.294284, -2.950189, 0.654665>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-4.254293, -1.245944, 1.556197>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<0.029664, -2.393158, 0.335538>, 0.272
	pigment { rgbt <0.400000005960464, 0.400000005960464, 0.400000005960464,0> }
}
sphere {
	<-0.622227, -0.364299, -2.012517>, 0.2432
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581,0> }
}
sphere {
	<-0.673779, -1.236077, -2.419907>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-2.321065, -3.247014, 0.240905>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-4.05502, -3.048759, -0.13346>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-3.552917, -3.661214, 1.453012>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-5.157072, -0.323684, -2.080738>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-5.032484, 0.817818, -0.738719>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<-5.40541, -0.895162, -0.421142>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.904345, -1.968492, -3.00088>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.45853, -0.921121, -2.745631>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.640345, -2.416988, -1.824519>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<4.336951, -2.506362, -1.042229>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<3.065234, -2.541784, 0.196345>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<4.455487, -1.431993, 0.364281>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<4.692695, -0.469255, -2.436925>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<4.537085, 0.741845, -1.15>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<3.453452, 0.788738, -2.563694>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<4.696633, 0.528244, 0.986588>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<4.31802, 1.670913, 2.286397>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<3.95127, 2.103295, 0.60632>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<3.310396, -1.294196, 2.191361>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.621977, -0.939404, 2.639461>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.971385, -0.132285, 3.491187>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<2.06319, 2.030449, 3.119241>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<0.638719, 1.372955, 2.265587>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<1.63082, 2.66104, 1.516873>, 0.176
	pigment { rgbt <0.75, 0.75, 0.75,0> }
}
sphere {
	<0.361959, -3.498898, 0.506357>, 0.2432
	pigment { rgbt <1, 0.0500000007450581, 0.0500000007450581,0> }
}
}
merge {
}
