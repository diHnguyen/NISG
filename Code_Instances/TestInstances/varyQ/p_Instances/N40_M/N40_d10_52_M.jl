global arcs = [1 7; 1 18; 1 25; 1 29; 1 37; 2 7; 2 29; 2 33; 2 38; 3 12; 3 19; 3 25; 3 40; 4 2; 4 14; 4 29; 5 2; 5 14; 5 22; 5 39; 6 28; 7 8; 7 24; 7 29; 7 30; 7 40; 8 6; 8 24; 8 35; 8 36; 9 10; 9 12; 9 15; 9 18; 9 40; 10 6; 10 12; 10 16; 10 30; 10 31; 10 33; 11 2; 11 10; 11 12; 11 15; 11 17; 11 32; 11 38; 12 18; 12 21; 12 27; 12 32; 12 38; 12 39; 13 28; 13 31; 14 6; 14 9; 14 28; 14 38; 15 21; 15 34; 15 36; 15 39; 16 26; 16 37; 17 11; 17 22; 17 28; 17 35; 17 39; 18 8; 18 12; 19 6; 19 12; 19 16; 19 17; 19 23; 19 25; 19 32; 19 37; 19 40; 20 13; 20 34; 21 11; 21 19; 21 28; 21 39; 22 18; 22 24; 22 38; 23 16; 23 21; 23 38; 24 3; 24 22; 24 27; 24 28; 24 35; 25 5; 25 10; 25 17; 25 18; 25 27; 25 36; 26 9; 26 18; 26 21; 26 36; 27 2; 27 4; 27 18; 27 33; 27 37; 28 6; 28 13; 28 29; 28 31; 29 3; 29 6; 29 13; 29 16; 29 39; 30 8; 30 37; 31 10; 31 19; 32 13; 32 22; 32 26; 32 34; 32 35; 33 36; 34 19; 34 26; 35 3; 35 5; 35 15; 35 20; 35 24; 35 31; 36 6; 36 23; 36 26; 37 3; 37 19; 37 32; 38 2; 38 19; 38 26; 38 36; 38 39; 39 19; 39 33; 39 40]
global d_x = [5.0, 8.0, 7.0, 1.0, 3.0, 4.0, 3.0, 1.0, 7.0, 3.0, 3.0, 3.0, 8.0, 5.0, 3.0, 10.0, 8.0, 5.0, 9.0, 4.0, 7.0, 2.0, 4.0, 7.0, 8.0, 3.0, 9.0, 6.0, 2.0, 3.0, 3.0, 7.0, 7.0, 9.0, 3.0, 5.0, 9.0, 4.0, 8.0, 9.0, 4.0, 1.0, 3.0, 5.0, 3.0, 9.0, 3.0, 1.0, 4.0, 2.0, 6.0, 2.0, 1.0, 1.0, 10.0, 8.0, 3.0, 6.0, 10.0, 2.0, 10.0, 8.0, 8.0, 10.0, 8.0, 3.0, 2.0, 7.0, 5.0, 5.0, 1.0, 8.0, 4.0, 10.0, 10.0, 4.0, 5.0, 10.0, 2.0, 3.0, 9.0, 6.0, 3.0, 4.0, 9.0, 2.0, 1.0, 7.0, 8.0, 6.0, 9.0, 1.0, 7.0, 6.0, 10.0, 1.0, 7.0, 10.0, 3.0, 10.0, 2.0, 4.0, 5.0, 7.0, 6.0, 2.0, 1.0, 1.0, 7.0, 5.0, 8.0, 2.0, 3.0, 9.0, 5.0, 10.0, 1.0, 7.0, 7.0, 8.0, 9.0, 7.0, 4.0, 3.0, 6.0, 1.0, 1.0, 5.0, 5.0, 1.0, 2.0, 2.0, 5.0, 1.0, 8.0, 10.0, 2.0, 8.0, 5.0, 8.0, 8.0, 4.0, 3.0, 3.0, 6.0, 1.0, 2.0, 2.0, 6.0, 5.0, 5.0, 9.0, 6.0, 1.0, 4.0]
global b_x = 5
global d_y = [6.0, 8.0, 2.0, 10.0, 10.0, 8.0, 5.0, 10.0, 4.0, 7.0, 2.0, 7.0, 4.0, 8.0, 3.0, 1.0, 6.0, 1.0, 5.0, 8.0, 10.0, 2.0, 7.0, 9.0, 1.0, 2.0, 6.0, 4.0, 8.0, 2.0, 6.0, 1.0, 1.0, 10.0, 6.0, 7.0, 9.0, 4.0, 4.0, 10.0, 4.0, 3.0, 4.0, 7.0, 3.0, 4.0, 1.0, 6.0, 8.0, 5.0, 8.0, 4.0, 5.0, 2.0, 8.0, 1.0, 1.0, 6.0, 7.0, 9.0, 9.0, 7.0, 9.0, 6.0, 10.0, 7.0, 4.0, 1.0, 10.0, 6.0, 1.0, 1.0, 4.0, 8.0, 2.0, 9.0, 6.0, 5.0, 10.0, 8.0, 5.0, 6.0, 8.0, 2.0, 7.0, 7.0, 10.0, 6.0, 1.0, 3.0, 6.0, 4.0, 4.0, 5.0, 10.0, 9.0, 5.0, 10.0, 4.0, 7.0, 8.0, 7.0, 9.0, 5.0, 9.0, 9.0, 9.0, 9.0, 2.0, 9.0, 2.0, 3.0, 4.0, 3.0, 2.0, 6.0, 10.0, 1.0, 7.0, 8.0, 10.0, 3.0, 3.0, 7.0, 7.0, 9.0, 5.0, 1.0, 10.0, 2.0, 3.0, 5.0, 8.0, 2.0, 2.0, 2.0, 8.0, 7.0, 10.0, 5.0, 2.0, 4.0, 4.0, 3.0, 8.0, 5.0, 4.0, 10.0, 6.0, 3.0, 10.0, 9.0, 4.0, 4.0, 3.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.645, 0.701, 0.621, 0.719, 0.669, 0.776, 0.704, 0.706, 0.613, 0.613, 0.79, 0.609, 0.639, 0.728, 0.755, 0.699, 0.781, 0.654, 0.733, 0.671, 0.716, 0.626, 0.745, 0.688, 0.628, 0.75, 0.786, 0.601, 0.666, 0.706, 0.623, 0.728, 0.705, 0.744, 0.754, 0.66, 0.735, 0.766, 0.719, 0.77, 0.695, 0.667, 0.605, 0.611, 0.733, 0.797, 0.686, 0.618, 0.601, 0.737, 0.701, 0.64, 0.742, 0.695, 0.763, 0.638, 0.748, 0.673, 0.723, 0.728, 0.755, 0.69, 0.7, 0.625, 0.71, 0.79, 0.696, 0.689, 0.799, 0.798, 0.724, 0.775, 0.733, 0.793, 0.709, 0.634, 0.705, 0.712, 0.701, 0.704, 0.76, 0.792, 0.709, 0.699, 0.718, 0.732, 0.63, 0.758, 0.681, 0.642, 0.779, 0.731, 0.638, 0.778, 0.725, 0.723, 0.785, 0.759, 0.67, 0.677, 0.792, 0.642, 0.679, 0.716, 0.714, 0.608, 0.786, 0.747, 0.726, 0.74, 0.611, 0.715, 0.675, 0.796, 0.665, 0.731, 0.616, 0.76, 0.795, 0.788, 0.747, 0.791, 0.643, 0.741, 0.622, 0.663, 0.618, 0.753, 0.618, 0.643, 0.736, 0.681, 0.645, 0.794, 0.688, 0.793, 0.723, 0.701, 0.624, 0.612, 0.777, 0.765, 0.653, 0.746, 0.79, 0.636, 0.633, 0.657, 0.721, 0.787, 0.783, 0.609, 0.672, 0.774, 0.694]
global origin = 1
global destination = 40