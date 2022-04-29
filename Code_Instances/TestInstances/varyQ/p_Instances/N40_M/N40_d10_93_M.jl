global arcs = [1 21; 1 22; 1 27; 1 38; 1 40; 2 23; 3 22; 3 28; 3 29; 4 13; 4 26; 5 8; 5 13; 5 14; 5 21; 5 26; 6 2; 6 20; 7 3; 7 12; 7 24; 7 25; 7 30; 8 16; 8 26; 8 29; 8 33; 8 37; 8 40; 9 23; 10 5; 10 16; 10 27; 10 35; 11 4; 11 6; 11 14; 12 8; 12 20; 12 22; 12 26; 13 9; 13 12; 13 14; 13 15; 13 18; 13 24; 13 28; 13 31; 13 37; 13 39; 14 33; 15 3; 15 25; 15 39; 15 40; 16 7; 16 11; 16 20; 16 28; 16 40; 17 7; 17 35; 17 36; 18 7; 18 30; 18 35; 19 10; 19 15; 19 27; 19 28; 19 37; 20 15; 20 18; 20 28; 21 3; 21 12; 21 22; 22 9; 22 38; 23 7; 23 27; 23 34; 24 7; 24 13; 24 15; 24 27; 25 13; 25 16; 25 23; 25 31; 25 38; 26 3; 26 12; 26 15; 26 16; 26 17; 26 23; 26 24; 27 17; 27 18; 27 35; 28 2; 28 3; 28 30; 28 40; 29 7; 29 9; 29 15; 29 24; 29 25; 29 33; 29 35; 29 37; 30 8; 30 19; 30 22; 31 9; 31 12; 31 24; 31 29; 31 30; 31 32; 31 39; 32 8; 32 13; 32 20; 32 28; 32 40; 33 7; 33 19; 34 20; 34 22; 34 28; 35 4; 36 8; 36 17; 36 28; 37 4; 37 14; 37 16; 37 39; 37 40; 38 10; 38 13; 38 14; 38 19; 38 22; 38 25; 38 30; 39 10; 39 21; 39 22; 39 27; 39 38]
global d_x = [6.0, 10.0, 8.0, 8.0, 8.0, 10.0, 8.0, 6.0, 2.0, 9.0, 2.0, 8.0, 6.0, 6.0, 3.0, 2.0, 4.0, 3.0, 8.0, 8.0, 3.0, 4.0, 1.0, 1.0, 3.0, 4.0, 6.0, 1.0, 5.0, 8.0, 9.0, 2.0, 2.0, 9.0, 5.0, 10.0, 7.0, 2.0, 8.0, 6.0, 6.0, 6.0, 6.0, 8.0, 8.0, 8.0, 8.0, 4.0, 6.0, 4.0, 5.0, 7.0, 1.0, 10.0, 5.0, 5.0, 2.0, 5.0, 6.0, 2.0, 8.0, 5.0, 10.0, 4.0, 5.0, 4.0, 4.0, 6.0, 2.0, 10.0, 4.0, 4.0, 2.0, 2.0, 9.0, 8.0, 3.0, 7.0, 7.0, 2.0, 3.0, 5.0, 6.0, 9.0, 4.0, 8.0, 5.0, 6.0, 2.0, 7.0, 8.0, 1.0, 10.0, 4.0, 2.0, 2.0, 7.0, 1.0, 1.0, 3.0, 4.0, 2.0, 8.0, 1.0, 7.0, 8.0, 8.0, 6.0, 7.0, 4.0, 6.0, 6.0, 8.0, 10.0, 2.0, 3.0, 4.0, 10.0, 5.0, 10.0, 6.0, 5.0, 2.0, 9.0, 5.0, 5.0, 7.0, 6.0, 6.0, 1.0, 1.0, 3.0, 5.0, 10.0, 3.0, 9.0, 4.0, 1.0, 5.0, 8.0, 9.0, 2.0, 2.0, 4.0, 7.0, 5.0, 7.0, 2.0, 8.0, 2.0, 3.0, 1.0, 5.0, 5.0, 3.0]
global b_x = 5
global d_y = [6.0, 2.0, 10.0, 9.0, 9.0, 10.0, 4.0, 8.0, 3.0, 8.0, 8.0, 8.0, 3.0, 10.0, 3.0, 3.0, 10.0, 7.0, 2.0, 1.0, 1.0, 4.0, 4.0, 8.0, 3.0, 6.0, 5.0, 3.0, 9.0, 2.0, 2.0, 2.0, 3.0, 8.0, 9.0, 3.0, 7.0, 10.0, 1.0, 5.0, 3.0, 6.0, 8.0, 8.0, 2.0, 4.0, 1.0, 1.0, 5.0, 1.0, 6.0, 4.0, 5.0, 7.0, 5.0, 8.0, 2.0, 5.0, 4.0, 6.0, 10.0, 3.0, 7.0, 5.0, 5.0, 10.0, 6.0, 8.0, 2.0, 1.0, 3.0, 4.0, 8.0, 6.0, 2.0, 10.0, 7.0, 2.0, 2.0, 10.0, 10.0, 4.0, 8.0, 2.0, 8.0, 3.0, 6.0, 6.0, 7.0, 10.0, 8.0, 7.0, 4.0, 2.0, 6.0, 7.0, 8.0, 5.0, 9.0, 2.0, 4.0, 1.0, 6.0, 1.0, 3.0, 1.0, 4.0, 1.0, 4.0, 3.0, 2.0, 10.0, 1.0, 2.0, 3.0, 9.0, 8.0, 6.0, 1.0, 7.0, 2.0, 5.0, 1.0, 1.0, 5.0, 10.0, 2.0, 10.0, 3.0, 6.0, 5.0, 7.0, 3.0, 7.0, 1.0, 4.0, 8.0, 8.0, 1.0, 1.0, 5.0, 6.0, 2.0, 6.0, 5.0, 7.0, 10.0, 4.0, 1.0, 4.0, 9.0, 3.0, 8.0, 4.0, 2.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.797, 0.719, 0.707, 0.648, 0.767, 0.797, 0.63, 0.668, 0.744, 0.701, 0.734, 0.643, 0.616, 0.609, 0.774, 0.76, 0.651, 0.657, 0.636, 0.718, 0.766, 0.658, 0.68, 0.746, 0.778, 0.748, 0.613, 0.707, 0.758, 0.652, 0.702, 0.713, 0.691, 0.615, 0.708, 0.604, 0.772, 0.764, 0.665, 0.717, 0.729, 0.639, 0.667, 0.725, 0.629, 0.628, 0.72, 0.614, 0.713, 0.714, 0.761, 0.643, 0.611, 0.787, 0.757, 0.623, 0.722, 0.728, 0.767, 0.653, 0.686, 0.756, 0.676, 0.609, 0.767, 0.718, 0.681, 0.759, 0.602, 0.755, 0.762, 0.707, 0.693, 0.788, 0.791, 0.661, 0.775, 0.745, 0.732, 0.621, 0.623, 0.648, 0.613, 0.648, 0.756, 0.683, 0.653, 0.783, 0.738, 0.717, 0.615, 0.681, 0.798, 0.787, 0.787, 0.647, 0.788, 0.759, 0.692, 0.788, 0.707, 0.704, 0.604, 0.69, 0.668, 0.755, 0.625, 0.769, 0.733, 0.758, 0.673, 0.635, 0.776, 0.769, 0.753, 0.741, 0.604, 0.748, 0.604, 0.686, 0.672, 0.734, 0.663, 0.704, 0.711, 0.667, 0.602, 0.75, 0.695, 0.61, 0.721, 0.689, 0.783, 0.788, 0.773, 0.772, 0.629, 0.723, 0.729, 0.684, 0.619, 0.703, 0.73, 0.654, 0.798, 0.615, 0.619, 0.648, 0.64, 0.778, 0.702, 0.788, 0.733, 0.75, 0.736]
global origin = 1
global destination = 40