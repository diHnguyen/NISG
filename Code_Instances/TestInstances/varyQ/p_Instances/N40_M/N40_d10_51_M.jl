global arcs = [1 31; 1 32; 1 39; 2 7; 2 20; 2 29; 2 34; 3 25; 3 38; 4 24; 4 39; 5 11; 5 17; 6 2; 6 3; 6 9; 6 12; 6 16; 6 18; 6 29; 7 10; 7 13; 7 30; 7 36; 8 9; 9 4; 9 17; 9 38; 10 31; 10 32; 10 38; 11 6; 11 18; 11 25; 11 37; 12 2; 12 11; 12 18; 12 38; 13 10; 13 17; 13 19; 14 17; 14 20; 14 25; 14 33; 14 34; 15 2; 15 7; 15 24; 15 31; 16 20; 16 27; 16 36; 17 21; 17 26; 17 34; 17 35; 17 39; 18 2; 18 5; 18 10; 18 19; 19 9; 19 29; 19 36; 20 2; 20 12; 20 19; 20 35; 21 8; 21 11; 21 13; 21 22; 21 30; 21 38; 21 40; 22 14; 22 30; 22 35; 22 38; 23 3; 23 22; 23 27; 24 10; 24 12; 24 23; 24 35; 25 32; 25 34; 26 17; 26 25; 26 28; 26 30; 27 23; 28 4; 28 17; 28 18; 29 15; 29 36; 29 40; 30 5; 30 6; 30 14; 30 27; 30 28; 30 29; 30 32; 30 35; 31 20; 32 8; 32 10; 32 19; 32 20; 33 5; 33 7; 33 9; 33 18; 33 23; 33 37; 33 38; 33 40; 34 15; 34 21; 34 29; 34 37; 35 3; 35 5; 35 8; 35 13; 35 26; 36 6; 36 18; 36 22; 36 23; 36 31; 36 32; 37 2; 37 34; 37 38; 38 14; 39 9; 39 21]
global d_x = [2.0, 4.0, 6.0, 10.0, 10.0, 10.0, 8.0, 10.0, 5.0, 4.0, 2.0, 1.0, 6.0, 2.0, 10.0, 1.0, 4.0, 7.0, 8.0, 2.0, 4.0, 6.0, 9.0, 8.0, 5.0, 2.0, 10.0, 3.0, 3.0, 3.0, 9.0, 2.0, 4.0, 5.0, 5.0, 8.0, 3.0, 4.0, 7.0, 4.0, 7.0, 2.0, 1.0, 9.0, 7.0, 2.0, 8.0, 2.0, 2.0, 6.0, 8.0, 7.0, 3.0, 1.0, 9.0, 4.0, 4.0, 6.0, 5.0, 10.0, 4.0, 7.0, 4.0, 7.0, 4.0, 2.0, 10.0, 4.0, 1.0, 8.0, 9.0, 5.0, 6.0, 5.0, 10.0, 8.0, 8.0, 8.0, 5.0, 7.0, 3.0, 5.0, 10.0, 8.0, 1.0, 3.0, 9.0, 7.0, 1.0, 4.0, 1.0, 1.0, 2.0, 8.0, 10.0, 5.0, 4.0, 2.0, 3.0, 6.0, 7.0, 2.0, 2.0, 5.0, 2.0, 9.0, 10.0, 5.0, 8.0, 9.0, 5.0, 2.0, 4.0, 4.0, 10.0, 2.0, 6.0, 7.0, 7.0, 1.0, 8.0, 10.0, 6.0, 4.0, 8.0, 7.0, 9.0, 6.0, 4.0, 4.0, 7.0, 9.0, 6.0, 8.0, 10.0, 6.0, 4.0, 1.0, 5.0, 2.0, 2.0, 4.0, 6.0]
global b_x = 5
global d_y = [7.0, 7.0, 9.0, 9.0, 5.0, 8.0, 2.0, 3.0, 1.0, 3.0, 7.0, 4.0, 9.0, 5.0, 1.0, 3.0, 1.0, 4.0, 4.0, 5.0, 9.0, 3.0, 3.0, 4.0, 1.0, 4.0, 7.0, 6.0, 9.0, 5.0, 5.0, 8.0, 10.0, 9.0, 2.0, 4.0, 6.0, 7.0, 10.0, 9.0, 4.0, 10.0, 9.0, 3.0, 5.0, 3.0, 5.0, 1.0, 2.0, 3.0, 7.0, 4.0, 10.0, 1.0, 9.0, 7.0, 10.0, 7.0, 8.0, 5.0, 10.0, 1.0, 7.0, 4.0, 6.0, 3.0, 9.0, 4.0, 10.0, 3.0, 5.0, 1.0, 3.0, 9.0, 9.0, 6.0, 1.0, 2.0, 2.0, 7.0, 5.0, 6.0, 8.0, 8.0, 1.0, 2.0, 7.0, 5.0, 9.0, 4.0, 9.0, 5.0, 8.0, 2.0, 9.0, 9.0, 2.0, 4.0, 6.0, 8.0, 5.0, 9.0, 6.0, 5.0, 5.0, 7.0, 5.0, 7.0, 7.0, 9.0, 8.0, 4.0, 4.0, 1.0, 4.0, 6.0, 8.0, 3.0, 7.0, 7.0, 9.0, 4.0, 9.0, 3.0, 9.0, 8.0, 9.0, 10.0, 8.0, 5.0, 6.0, 1.0, 6.0, 2.0, 5.0, 4.0, 4.0, 9.0, 3.0, 4.0, 10.0, 6.0, 3.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.639, 0.626, 0.786, 0.786, 0.775, 0.679, 0.696, 0.661, 0.794, 0.715, 0.647, 0.609, 0.734, 0.769, 0.647, 0.726, 0.777, 0.749, 0.761, 0.684, 0.615, 0.651, 0.791, 0.705, 0.679, 0.602, 0.677, 0.788, 0.636, 0.781, 0.71, 0.708, 0.616, 0.79, 0.68, 0.727, 0.707, 0.775, 0.709, 0.753, 0.791, 0.692, 0.62, 0.734, 0.735, 0.721, 0.666, 0.618, 0.742, 0.63, 0.649, 0.677, 0.722, 0.63, 0.682, 0.686, 0.655, 0.745, 0.7, 0.787, 0.688, 0.754, 0.792, 0.744, 0.657, 0.697, 0.788, 0.612, 0.759, 0.745, 0.786, 0.677, 0.72, 0.793, 0.606, 0.65, 0.724, 0.752, 0.603, 0.768, 0.784, 0.626, 0.738, 0.681, 0.646, 0.746, 0.755, 0.744, 0.655, 0.654, 0.757, 0.741, 0.71, 0.612, 0.789, 0.723, 0.753, 0.665, 0.629, 0.777, 0.64, 0.607, 0.76, 0.671, 0.793, 0.679, 0.719, 0.718, 0.773, 0.764, 0.798, 0.608, 0.638, 0.702, 0.776, 0.698, 0.683, 0.684, 0.777, 0.626, 0.661, 0.702, 0.799, 0.655, 0.615, 0.631, 0.602, 0.726, 0.603, 0.651, 0.748, 0.655, 0.737, 0.741, 0.685, 0.631, 0.762, 0.676, 0.672, 0.741, 0.632, 0.754, 0.672]
global origin = 1
global destination = 40