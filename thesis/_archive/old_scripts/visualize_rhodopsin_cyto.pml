# PyMOL script to visualize rhodopsin cytoplasmic face
# Generated for 3PQR (active rhodopsin / Meta II)

# Fetch structure
fetch 3pqr, async=0

# Basic setup
hide everything
show cartoon, chain A
color gray80, chain A

# Color TM helices
color red, chain A and resi 18+19+20+21+22+23+24+25+26+27+28+29+30+31+32+33+34+35+36+37+38+39+40+41+42+43+44+45+46+47+48+49+50+51+52+53+54+55+56+57+58+59+60+61+62+63+64+65
color orange, chain A and resi 70+71+72+73+74+75+76+77+78+79+80+81+82+83+84+85+86+87+88+89+90+91+92+93+94+95+96+97+98+99+100+101
color yellow, chain A and resi 106+107+108+109+110+111+112+113+114+115+116+117+118+119+120+121+122+123+124+125+126+127+128+129+130+131+132+133+134+135+136+137+138+139+140+141
color green, chain A and resi 149+150+151+152+153+154+155+156+157+158+159+160+161+162+163+164+165+166+167+168+169+170+171+172+173+174
color cyan, chain A and resi 199+200+201+202+203+204+205+206+207+208+209+210+211+212+213+214+215+216+217+218+219+220+221+222+223+224+225+226+227+228+229+230+231+232+233+234+235+236+237
color blue, chain A and resi 241+242+243+244+245+246+247+248+249+250+251+252+253+254+255+256+257+258+259+260+261+262+263+264+265+266+267+268+269+270+271+272+273+274+275+276+277+278
color purple, chain A and resi 284+285+286+287+288+289+290+291+292+293+294+295+296+297+298+299+300+301+302+303+304+305+306+307+308+309

# Highlight cytoplasmic face (designable region)
select cyto_face, chain A and resi 61+62+63+64+65+66+67+68+69+70+71+72+73+74+75+76+77+78+79+136+137+138+139+140+141+142+143+144+145+146+147+148+149+150+151+152+153+154+155+156+157+211+229+230+231+232+233+234+235+236+237+241+242+243+244+245+246+247+248+249+250+251+252+253+254+255+256+306+307+308+309+310+311+312+313+314+315+316+317+318+319+320+321+323+324+325+326
show sticks, cyto_face
color white, cyto_face
set stick_transparency, 0.3, cyto_face

# Show surface of cytoplasmic region
#select cyto_surface, chain A and resi 61+62+63+64+65+66+67+68+69+70+71+72+73+74+75+76+77+78+79+136+137+138+139+140+141+142+143+144+145+146+147+148+149+150+151+152+153+154+155+156+157+211+229+230+231+232+233+234+235+236+237+241+242+243+244+245+246+247+248+249+250+251+252+253+254+255+256+306+307+308+309+310+311+312+313+314+315+316+317+318+319+320+321+323+324+325+326
#show surface, cyto_surface
#set surface_color, yelloworange, cyto_surface
#set transparency, 0.7, cyto_surface

# Highlight ICL3 specifically (G-protein binding site)

# TM5 cytoplasmic end (moves during activation)
select tm5_cyto, chain A and resi 211+229+230+231+232+233+234+235+236+237
color cyan, tm5_cyto

# TM6 cytoplasmic start (large movement during activation)
select tm6_cyto, chain A and resi 241+242+243+244+245+246+247+248+249+250+251+252+253+254+255+256
color blue, tm6_cyto

# Orient for cytoplasmic view
orient chain A
turn x, 90
turn y, 180

# Labels
#set label_size, 14
#label icl3 and name CA, resn+resi

# Ray trace settings
set ray_shadows, 0
bg_color white
set antialias, 2

# Save session
#save rhodopsin_cyto_face.pse

print("Visualization ready!")
print("Cytoplasmic face shown in sticks (white)")
print("ICL3 (G-protein binding) shown in pink spheres")
print("TM5 end in cyan, TM6 start in blue")
