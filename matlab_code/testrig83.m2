KK = ZZ / 30097;
R = KK[x1,x2,x3,x4,x5,x6,x7];
eq1 = 10000*x1*x3*x6-16777999432*x1*x3+2500*x1*x5^2-10000*x2^2*x6+16777999432*x2^2-5000*x2*x4*x5+2500*x3*x4^2
eq2 = 345696765*x1*x3-345696765*x2^2-5000*x2*x5+5000*x3*x4+2500*x3
eq3 = -4945732904*x1*x3+5000*x1*x5+2500*x1+4945732904*x2^2-5000*x2*x4
eq4 = -1.262056216723612e+18*x1^2+2.036721962299979e+19*x1*x2-3.183982650795666e+19*x1*x3-4119687717463551*x1*x4+1.782087531260572e+16*x1*x5-5141125275849*x1*x6-224682550878*x1*x7+1.130031807569803e+17*x1-8.217217855960873e+19*x2^2+2.569175329326006e+20*x2*x3+3.324201545378951e+16*x2*x4-1.437977471036156e+17*x2*x5+35662698866020*x2*x6+1812977424676*x2*x7-9.118296673621942e+17*x2-2.008180258975727e+20*x3^2-5.196683810627966e+16*x3*x4+2.247972676075953e+17*x3*x5-61845803309979*x3*x6-2834205538754*x3*x7+1.425452219351285e+18*x3-3361939560323*x4^2+29086042351737*x4*x5+226740497*x4*x6-366711852*x4*x7+184436243659341*x4-62909954544638*x5^2-786421028*x5*x6+1586315915*x5*x7-797831177214118*x5+2604504368*x6+10058919156*x7-2529546364888931
eq5 = 9.06615095161258e+21*x1^2-4.478197120649658e+22*x1*x2+1.565824493143516e+25*x1*x3+5.60058203913436e+18*x1*x4-7.307433231397776e+19*x1*x5-914021074503666*x1*x6+806906353873576*x1*x7-5.201828422374392e+20*x1-1.575282488530684e+25*x2^2+1.964986709433123e+23*x2*x3+5.879337510038424e+19*x2*x4-8.442334833502616e+19*x2*x5-2440330129809385*x2*x6+2526546549540044*x2*x7-1.530145611781533e+21*x2+2.523310409247726e+23*x3^2-4.21889317308376e+18*x3*x4-7.200643290506824e+19*x3*x5-1798691098181870*x3*x6+1779195047439589*x3*x7-1.15342569263055e+21*x3-56247106410332*x4^2+150745442311475*x4*x5+2578666256126695*x4+219163775891897*x5^2+3633982516450677*x5+991020727949858
eq6 = 7.260467162878771e+20*x1^2-2.440423709253829e+22*x1*x2+2.23705014111517e+25*x1*x3+4.669314426703469e+17*x1*x4-3.886922705424106e+19*x1*x5-123407608906862*x1*x6+64516328236401*x1*x7-1.779374528667242e+20*x1-2.213202159745439e+25*x2^2-8.768407083094219e+23*x2*x3+2.44835821758288e+19*x2*x4+5.61992259339168e+19*x2*x5+1616470390932357*x2*x6-1649931200239833*x2*x7+1.019444020865302e+21*x2+1.005037810041162e+24*x3^2+6.068898666244176e+19*x3*x4-2.852823302751776e+20*x3*x5-7098920566210090*x3*x6+7090784195255637*x3*x7-4.492453860471829e+21*x3-5556704169449*x4^2-141455677636804*x4*x5-1134597983384210*x4+459145822905129*x5^2+8805895388911655*x5+2.417317019194895e+16
eq7 = 1.26205621672361e+18*x1^2-2.036721962299978e+19*x1*x2+5.033235205165214e+19*x1*x3+1572456060283544*x1*x4-3702124199841808*x1*x5-1734515485436*x1*x6-1.280296818145402e+18*x1+6.367965301591329e+19*x2^2-2.569175329326007e+20*x2*x3-1.797224236429771e+16*x2*x4+3.967075910516086e+16*x2*x5+20579496894726*x2*x6+1.520212364246899e+19*x2+2.008180258975727e+20*x3^2+5.267418905873522e+16*x3*x4-1.133533259563798e+17*x3*x5-61148977484502*x3*x6-4.518343374937774e+19*x3+112341275439*x4^2-906488712338*x4*x5-139971356*x4*x6-101632175239650*x4+1417102769377*x5^2+799894887*x5*x6+587301803876419*x5+12663423525*x6+3897015137978983
eqs = {eq1,eq2,eq3,eq4,eq5,eq6,eq7};
I = ideal eqs;
gbTrace = 3
dim I, degree I
