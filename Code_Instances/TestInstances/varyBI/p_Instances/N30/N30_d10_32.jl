global arcs = [1 2; 1 4; 1 5; 1 14; 1 16; 1 19; 2 10; 2 14; 2 22; 2 27; 3 17; 4 20; 5 3; 5 12; 5 29; 6 4; 6 5; 7 2; 7 10; 7 19; 8 26; 9 16; 9 20; 9 22; 9 28; 10 3; 10 12; 10 18; 10 20; 11 6; 11 7; 11 9; 11 26; 11 27; 12 6; 12 7; 12 27; 13 4; 13 10; 13 17; 13 19; 13 24; 14 15; 15 7; 15 11; 15 16; 15 21; 16 14; 16 17; 17 4; 17 5; 17 26; 18 6; 18 7; 18 8; 18 15; 18 20; 18 30; 19 4; 19 8; 19 23; 20 3; 20 15; 21 4; 21 20; 22 3; 22 4; 22 8; 22 11; 22 24; 23 4; 23 20; 23 30; 24 13; 24 20; 24 21; 25 2; 25 15; 25 22; 25 29; 26 3; 27 9; 27 25; 28 2; 28 4; 28 17; 28 27; 29 7; 29 11; 29 14; 29 17; 29 26]
global d_x = [8.0, 5.0, 9.0, 8.0, 9.0, 6.0, 2.0, 5.0, 5.0, 1.0, 3.0, 8.0, 8.0, 6.0, 4.0, 10.0, 6.0, 10.0, 6.0, 8.0, 6.0, 1.0, 5.0, 2.0, 4.0, 7.0, 1.0, 9.0, 10.0, 1.0, 7.0, 10.0, 2.0, 1.0, 10.0, 4.0, 8.0, 8.0, 8.0, 9.0, 5.0, 1.0, 2.0, 3.0, 5.0, 1.0, 9.0, 6.0, 1.0, 7.0, 4.0, 5.0, 9.0, 2.0, 8.0, 2.0, 5.0, 6.0, 6.0, 3.0, 9.0, 5.0, 10.0, 4.0, 5.0, 3.0, 4.0, 7.0, 3.0, 4.0, 4.0, 3.0, 9.0, 3.0, 10.0, 8.0, 4.0, 9.0, 7.0, 8.0, 10.0, 10.0, 2.0, 9.0, 7.0, 6.0, 2.0, 8.0, 5.0, 4.0, 7.0, 4.0]
global b_x = 5
global d_y = [7.0, 10.0, 10.0, 3.0, 1.0, 7.0, 7.0, 3.0, 10.0, 7.0, 4.0, 3.0, 10.0, 3.0, 7.0, 10.0, 5.0, 1.0, 7.0, 2.0, 5.0, 5.0, 1.0, 6.0, 8.0, 8.0, 10.0, 6.0, 4.0, 3.0, 2.0, 3.0, 6.0, 2.0, 1.0, 8.0, 9.0, 5.0, 6.0, 4.0, 6.0, 4.0, 3.0, 10.0, 8.0, 8.0, 1.0, 3.0, 5.0, 3.0, 4.0, 2.0, 4.0, 4.0, 6.0, 10.0, 7.0, 3.0, 7.0, 9.0, 8.0, 5.0, 10.0, 4.0, 10.0, 9.0, 10.0, 4.0, 5.0, 3.0, 2.0, 10.0, 5.0, 3.0, 9.0, 3.0, 3.0, 8.0, 5.0, 4.0, 3.0, 6.0, 9.0, 8.0, 1.0, 3.0, 1.0, 6.0, 4.0, 9.0, 6.0, 2.0]
global b_y = 10
global p = [0.023, 0.015, 0.138, 0.461, 0.456, 0.696, 0.029, 0.635, 0.115, 0.732, 0.425, 0.308, 0.398, 0.142, 0.512, 0.411, 0.833, 0.323, 0.419, 0.801, 0.962, 0.897, 0.424, 0.512, 0.658, 0.567, 0.387, 0.328, 0.826, 0.587, 0.1, 0.733, 0.401, 0.016, 0.609, 0.592, 0.096, 0.125, 0.118, 0.157, 0.036, 0.015, 0.384, 0.67, 0.476, 0.288, 0.037, 0.646, 0.221, 0.081, 0.321, 0.178, 0.203, 0.764, 0.884, 0.266, 0.927, 0.324, 0.802, 0.412, 0.321, 0.608, 0.136, 0.523, 0.095, 0.655, 0.431, 0.711, 0.088, 0.02, 0.739, 0.595, 0.179, 0.636, 0.401, 0.463, 0.287, 0.62, 0.88, 0.49, 0.436, 0.787, 0.367, 0.779, 0.007, 0.199, 0.028, 0.943, 0.4, 0.776, 0.673, 0.242]
global q = [0.038, 0.339, 0.632, 0.691, 0.698, 0.772, 0.435, 0.767, 0.164, 0.968, 0.761, 0.969, 0.973, 0.633, 0.62, 0.475, 0.855, 0.864, 0.682, 0.942, 0.965, 0.994, 0.732, 0.914, 0.847, 0.679, 0.471, 0.876, 0.948, 0.751, 0.113, 0.827, 0.672, 0.596, 0.936, 0.667, 0.595, 0.157, 0.686, 0.722, 0.485, 0.454, 0.471, 0.837, 0.554, 0.577, 0.219, 0.785, 0.967, 0.767, 0.823, 0.561, 0.205, 0.895, 0.987, 0.472, 0.929, 0.575, 0.877, 0.519, 0.402, 0.888, 0.594, 0.893, 0.265, 0.889, 0.985, 0.888, 0.176, 0.52, 0.947, 0.747, 0.611, 0.989, 0.7, 0.618, 0.595, 0.675, 0.9, 0.498, 0.755, 0.902, 0.914, 0.789, 0.168, 0.66, 0.771, 0.997, 0.565, 0.927, 0.844, 0.735]
global origin = 1
global destination = 30
