global arcs = [1 6; 1 19; 2 5; 2 6; 2 9; 2 15; 2 17; 2 26; 2 37; 3 9; 3 34; 3 36; 3 39; 3 40; 4 9; 4 26; 4 31; 5 3; 5 10; 5 30; 5 32; 5 33; 6 11; 6 16; 6 23; 7 2; 7 13; 7 15; 7 29; 7 35; 8 11; 8 12; 8 15; 8 24; 8 27; 9 15; 9 16; 9 21; 9 30; 9 32; 9 33; 9 34; 10 5; 10 7; 10 14; 10 16; 10 20; 10 37; 11 4; 11 9; 11 10; 11 14; 11 17; 11 28; 11 30; 12 18; 12 38; 13 24; 13 31; 13 35; 14 2; 14 9; 14 13; 14 17; 14 20; 14 34; 15 5; 15 8; 15 26; 16 17; 16 36; 17 23; 18 9; 18 16; 18 21; 18 34; 18 40; 19 6; 19 25; 19 26; 19 28; 20 12; 20 37; 21 18; 21 26; 21 32; 22 11; 22 23; 23 12; 23 16; 23 22; 23 25; 23 31; 23 35; 24 10; 24 29; 24 35; 25 3; 25 13; 25 17; 25 28; 25 35; 26 15; 26 22; 27 7; 27 13; 27 24; 28 5; 28 30; 28 32; 28 33; 29 9; 29 34; 30 4; 30 19; 30 23; 30 24; 30 31; 30 33; 30 36; 31 9; 31 11; 31 19; 31 38; 32 18; 32 19; 32 20; 32 23; 33 7; 33 16; 33 17; 33 28; 33 35; 33 38; 34 9; 34 10; 34 17; 34 33; 34 39; 35 4; 35 5; 35 7; 35 18; 36 16; 36 38; 37 12; 37 20; 38 3; 38 10; 38 22; 38 36; 39 12; 39 26; 39 32; 39 33]
global d_x = [6.0, 10.0, 6.0, 2.0, 8.0, 1.0, 3.0, 3.0, 2.0, 5.0, 10.0, 4.0, 4.0, 9.0, 5.0, 9.0, 5.0, 1.0, 3.0, 2.0, 3.0, 1.0, 1.0, 8.0, 2.0, 1.0, 1.0, 4.0, 3.0, 4.0, 8.0, 9.0, 3.0, 10.0, 6.0, 5.0, 8.0, 6.0, 5.0, 5.0, 7.0, 10.0, 9.0, 6.0, 1.0, 9.0, 6.0, 4.0, 3.0, 10.0, 10.0, 10.0, 8.0, 6.0, 10.0, 6.0, 6.0, 8.0, 9.0, 3.0, 3.0, 2.0, 8.0, 5.0, 3.0, 5.0, 5.0, 10.0, 10.0, 1.0, 3.0, 9.0, 4.0, 10.0, 10.0, 8.0, 6.0, 1.0, 10.0, 2.0, 7.0, 1.0, 6.0, 10.0, 10.0, 9.0, 9.0, 7.0, 2.0, 8.0, 9.0, 4.0, 2.0, 2.0, 2.0, 5.0, 1.0, 8.0, 4.0, 9.0, 3.0, 6.0, 6.0, 3.0, 4.0, 3.0, 2.0, 5.0, 7.0, 10.0, 7.0, 5.0, 2.0, 3.0, 5.0, 1.0, 2.0, 8.0, 6.0, 6.0, 7.0, 7.0, 3.0, 8.0, 9.0, 2.0, 2.0, 9.0, 7.0, 9.0, 1.0, 5.0, 10.0, 2.0, 10.0, 4.0, 8.0, 9.0, 7.0, 5.0, 2.0, 8.0, 8.0, 1.0, 10.0, 9.0, 9.0, 10.0, 1.0, 1.0, 9.0, 5.0, 4.0, 4.0, 2.0]
global b_x = 5
global d_y = [5.0, 8.0, 4.0, 4.0, 9.0, 4.0, 5.0, 8.0, 5.0, 2.0, 9.0, 6.0, 9.0, 5.0, 2.0, 6.0, 6.0, 10.0, 6.0, 4.0, 10.0, 3.0, 5.0, 6.0, 10.0, 9.0, 5.0, 1.0, 4.0, 10.0, 9.0, 9.0, 4.0, 6.0, 2.0, 1.0, 10.0, 4.0, 6.0, 8.0, 7.0, 4.0, 1.0, 8.0, 1.0, 5.0, 8.0, 3.0, 9.0, 10.0, 3.0, 2.0, 10.0, 4.0, 6.0, 9.0, 6.0, 8.0, 8.0, 1.0, 6.0, 9.0, 1.0, 4.0, 10.0, 1.0, 3.0, 8.0, 8.0, 5.0, 5.0, 4.0, 3.0, 1.0, 10.0, 2.0, 8.0, 8.0, 10.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.0, 6.0, 4.0, 5.0, 10.0, 1.0, 4.0, 1.0, 10.0, 9.0, 2.0, 5.0, 1.0, 9.0, 4.0, 5.0, 2.0, 2.0, 6.0, 10.0, 5.0, 5.0, 9.0, 8.0, 3.0, 2.0, 2.0, 9.0, 7.0, 7.0, 4.0, 5.0, 9.0, 3.0, 5.0, 6.0, 4.0, 6.0, 2.0, 10.0, 7.0, 2.0, 8.0, 5.0, 9.0, 7.0, 3.0, 6.0, 2.0, 9.0, 5.0, 3.0, 10.0, 9.0, 2.0, 4.0, 3.0, 3.0, 1.0, 8.0, 5.0, 7.0, 8.0, 1.0, 3.0, 4.0, 8.0, 10.0, 3.0, 10.0, 2.0]
global b_y = 10
global p = [0.587, 0.23, 0.804, 0.261, 0.938, 0.889, 0.512, 0.631, 0.983, 0.78, 0.567, 0.518, 0.89, 0.284, 0.276, 0.46, 0.427, 0.335, 0.252, 0.699, 0.682, 0.596, 0.15, 0.27, 0.19, 0.25, 0.192, 0.703, 0.093, 0.588, 0.833, 0.885, 0.186, 0.538, 0.383, 0.383, 0.116, 0.968, 0.472, 0.264, 0.534, 0.347, 0.913, 0.942, 0.123, 0.667, 0.8, 0.779, 0.921, 0.82, 0.108, 0.002, 0.572, 0.91, 0.78, 0.945, 0.994, 0.499, 0.106, 0.181, 0.186, 0.925, 0.4, 0.577, 0.734, 0.549, 0.658, 0.218, 0.294, 0.237, 0.313, 0.767, 0.436, 0.791, 0.168, 0.291, 0.107, 0.767, 0.134, 0.908, 0.217, 0.292, 0.115, 0.18, 0.618, 0.722, 0.684, 0.802, 0.92, 0.07, 0.9, 0.491, 0.34, 0.234, 0.606, 0.625, 0.334, 0.019, 0.605, 0.925, 0.998, 0.279, 0.186, 0.565, 0.327, 0.753, 0.43, 0.125, 0.677, 0.227, 0.852, 0.514, 0.795, 0.802, 0.026, 0.513, 0.983, 0.801, 0.571, 0.81, 0.948, 0.861, 0.078, 0.67, 0.243, 0.495, 0.172, 0.775, 0.124, 0.513, 0.869, 0.783, 0.264, 0.553, 0.772, 0.358, 0.367, 0.653, 0.178, 0.092, 0.051, 0.581, 0.109, 0.882, 0.786, 0.732, 0.109, 0.86, 0.52, 0.783, 0.432, 0.513, 0.675, 0.402, 0.98]
global q = [0.94, 0.719, 0.872, 0.84, 0.994, 0.979, 0.997, 0.634, 0.984, 0.811, 0.647, 0.854, 0.976, 0.436, 0.279, 0.536, 0.583, 0.393, 0.728, 0.726, 0.809, 0.905, 0.31, 0.51, 0.982, 0.957, 0.547, 0.81, 0.378, 0.646, 0.972, 0.978, 0.89, 0.757, 0.718, 0.475, 0.163, 0.987, 0.844, 0.651, 0.737, 0.462, 0.996, 0.972, 0.977, 0.775, 0.838, 0.811, 0.965, 0.976, 0.365, 0.493, 0.832, 0.982, 0.999, 0.976, 0.999, 0.541, 0.756, 0.559, 0.323, 0.976, 0.626, 0.927, 0.947, 0.804, 0.66, 0.34, 0.417, 0.461, 0.833, 0.945, 0.935, 0.95, 0.638, 0.794, 0.257, 0.958, 0.338, 0.991, 0.886, 0.448, 0.729, 0.42, 0.705, 0.852, 0.754, 0.806, 0.988, 0.793, 0.957, 0.606, 0.458, 0.306, 0.842, 0.716, 0.603, 0.75, 0.896, 0.937, 0.998, 0.592, 0.984, 0.674, 0.386, 0.77, 0.735, 0.136, 0.909, 0.822, 0.956, 0.623, 0.857, 0.893, 0.653, 0.627, 0.995, 0.875, 0.576, 0.954, 0.963, 0.997, 0.342, 0.756, 0.505, 0.725, 0.378, 0.825, 0.531, 0.673, 0.905, 0.987, 0.709, 0.95, 0.966, 0.768, 0.46, 0.703, 0.305, 0.81, 0.094, 0.724, 0.756, 0.992, 0.969, 0.882, 0.576, 0.957, 0.798, 0.859, 0.771, 0.564, 0.783, 0.877, 0.993]
global origin = 1
global destination = 40