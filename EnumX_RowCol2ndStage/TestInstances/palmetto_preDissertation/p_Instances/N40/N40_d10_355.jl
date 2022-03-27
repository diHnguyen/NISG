global arcs = [1 2; 1 10; 2 10; 2 13; 2 19; 2 23; 2 33; 3 14; 3 36; 3 38; 4 8; 4 23; 4 30; 5 2; 5 14; 5 32; 6 12; 6 34; 7 10; 7 18; 7 24; 7 27; 8 15; 8 23; 8 26; 8 32; 8 37; 9 15; 9 28; 9 29; 9 32; 10 18; 10 20; 10 32; 10 33; 11 8; 12 8; 12 30; 12 35; 13 5; 13 17; 13 34; 13 38; 14 5; 14 21; 14 37; 15 24; 16 13; 16 27; 16 29; 16 30; 16 39; 17 9; 17 14; 17 30; 18 2; 18 21; 18 25; 18 29; 18 33; 18 39; 19 4; 19 5; 19 16; 19 23; 19 29; 19 32; 19 33; 20 12; 20 15; 20 18; 20 29; 21 3; 21 22; 21 35; 22 11; 22 24; 23 6; 23 22; 23 34; 24 13; 24 14; 24 25; 25 10; 25 33; 25 40; 26 7; 26 23; 27 9; 27 13; 27 15; 27 22; 27 29; 27 32; 27 33; 28 21; 28 23; 28 25; 28 32; 29 5; 29 7; 29 13; 29 15; 29 25; 29 30; 29 36; 29 40; 30 6; 30 11; 30 13; 30 18; 30 40; 31 6; 31 13; 31 19; 31 21; 31 26; 32 13; 32 15; 32 16; 32 25; 33 3; 33 23; 33 36; 33 38; 34 4; 34 15; 34 16; 34 21; 34 33; 35 19; 35 28; 35 32; 36 7; 36 8; 36 20; 36 22; 36 27; 36 30; 36 39; 37 4; 37 8; 37 12; 37 31; 38 16; 38 20; 38 26; 38 32; 39 18; 39 23]
global d_x = [4.0, 2.0, 10.0, 2.0, 10.0, 10.0, 7.0, 2.0, 10.0, 2.0, 8.0, 7.0, 6.0, 1.0, 1.0, 7.0, 6.0, 9.0, 10.0, 3.0, 6.0, 10.0, 6.0, 1.0, 5.0, 8.0, 3.0, 1.0, 1.0, 9.0, 7.0, 4.0, 9.0, 3.0, 2.0, 10.0, 1.0, 1.0, 2.0, 7.0, 5.0, 8.0, 6.0, 9.0, 7.0, 5.0, 4.0, 3.0, 7.0, 8.0, 6.0, 1.0, 2.0, 5.0, 3.0, 4.0, 7.0, 9.0, 5.0, 6.0, 9.0, 3.0, 1.0, 6.0, 4.0, 7.0, 7.0, 9.0, 5.0, 10.0, 7.0, 7.0, 5.0, 3.0, 9.0, 10.0, 10.0, 6.0, 9.0, 1.0, 1.0, 5.0, 5.0, 1.0, 1.0, 5.0, 7.0, 10.0, 7.0, 6.0, 1.0, 4.0, 9.0, 9.0, 5.0, 7.0, 5.0, 3.0, 5.0, 1.0, 3.0, 7.0, 2.0, 6.0, 8.0, 6.0, 3.0, 4.0, 6.0, 1.0, 6.0, 6.0, 10.0, 5.0, 4.0, 9.0, 5.0, 7.0, 5.0, 2.0, 4.0, 6.0, 9.0, 4.0, 8.0, 3.0, 3.0, 3.0, 6.0, 3.0, 7.0, 5.0, 2.0, 6.0, 2.0, 5.0, 1.0, 5.0, 7.0, 4.0, 6.0, 2.0, 10.0, 1.0, 6.0, 6.0, 9.0, 9.0, 6.0, 3.0]
global b_x = 5
global d_y = [4.0, 4.0, 1.0, 2.0, 6.0, 4.0, 7.0, 1.0, 5.0, 1.0, 5.0, 5.0, 5.0, 5.0, 2.0, 8.0, 9.0, 5.0, 6.0, 1.0, 6.0, 4.0, 2.0, 2.0, 2.0, 8.0, 3.0, 6.0, 6.0, 3.0, 1.0, 10.0, 4.0, 9.0, 1.0, 2.0, 1.0, 7.0, 3.0, 1.0, 3.0, 7.0, 3.0, 9.0, 3.0, 10.0, 5.0, 4.0, 3.0, 3.0, 4.0, 1.0, 7.0, 8.0, 6.0, 10.0, 6.0, 2.0, 2.0, 9.0, 2.0, 3.0, 10.0, 10.0, 3.0, 4.0, 3.0, 5.0, 9.0, 5.0, 3.0, 7.0, 6.0, 2.0, 1.0, 5.0, 3.0, 7.0, 4.0, 4.0, 1.0, 6.0, 1.0, 9.0, 6.0, 9.0, 5.0, 2.0, 1.0, 5.0, 9.0, 9.0, 4.0, 9.0, 7.0, 10.0, 4.0, 6.0, 6.0, 10.0, 10.0, 7.0, 3.0, 1.0, 1.0, 3.0, 2.0, 10.0, 7.0, 4.0, 8.0, 7.0, 5.0, 4.0, 1.0, 3.0, 2.0, 7.0, 6.0, 10.0, 4.0, 7.0, 9.0, 3.0, 9.0, 6.0, 8.0, 4.0, 4.0, 2.0, 4.0, 3.0, 1.0, 10.0, 10.0, 5.0, 1.0, 3.0, 10.0, 3.0, 7.0, 5.0, 1.0, 6.0, 8.0, 9.0, 2.0, 9.0, 6.0, 8.0]
global b_y = 10
global p = [0.348, 0.555, 0.706, 0.81, 0.664, 0.786, 0.836, 0.408, 0.148, 0.994, 0.993, 0.358, 0.992, 0.03, 0.807, 0.082, 0.753, 0.175, 0.236, 0.947, 0.755, 0.992, 0.873, 0.148, 0.79, 0.313, 0.496, 0.017, 0.779, 0.766, 0.372, 0.754, 0.65, 0.545, 0.206, 0.838, 0.025, 0.97, 0.946, 0.421, 0.208, 0.463, 0.68, 0.485, 0.972, 0.14, 0.337, 0.165, 0.897, 0.11, 0.364, 0.02, 0.878, 0.206, 0.828, 0.469, 0.945, 0.207, 0.098, 0.77, 0.291, 0.847, 0.31, 0.945, 0.098, 0.488, 0.741, 0.689, 0.642, 0.259, 0.351, 0.207, 0.935, 0.211, 0.649, 0.339, 0.171, 0.157, 0.685, 0.684, 0.763, 0.502, 0.507, 0.227, 0.59, 0.35, 0.567, 0.903, 0.934, 0.879, 0.878, 0.814, 0.455, 0.236, 0.713, 0.371, 0.719, 0.77, 0.428, 0.721, 0.297, 0.256, 0.47, 0.39, 0.688, 0.563, 0.713, 0.353, 0.377, 0.604, 0.174, 0.389, 0.818, 0.216, 0.558, 0.758, 0.899, 0.989, 0.51, 0.76, 0.768, 0.207, 0.69, 0.259, 0.505, 0.225, 0.415, 0.102, 0.278, 0.198, 0.629, 0.072, 0.704, 0.034, 0.564, 0.319, 0.476, 0.398, 0.938, 0.232, 0.069, 0.753, 0.717, 0.267, 0.468, 0.635, 0.717, 0.359, 0.908, 0.414]
global q = [0.38, 0.581, 0.966, 0.927, 0.876, 0.834, 0.847, 0.896, 0.552, 0.996, 0.997, 0.47, 0.996, 0.196, 0.914, 0.752, 0.962, 0.544, 0.875, 0.994, 0.901, 0.993, 0.891, 0.182, 0.901, 0.431, 0.879, 0.64, 0.786, 0.833, 0.993, 0.919, 0.938, 0.991, 0.839, 0.978, 0.891, 0.99, 0.952, 0.855, 0.347, 0.661, 0.726, 0.631, 0.989, 0.373, 0.924, 0.358, 0.981, 0.398, 0.937, 0.956, 0.902, 0.58, 0.93, 0.63, 0.961, 0.256, 0.648, 0.969, 0.515, 0.968, 0.727, 0.96, 0.243, 0.7, 0.905, 0.714, 0.99, 0.909, 0.514, 0.327, 0.981, 0.674, 0.841, 0.748, 0.805, 0.991, 0.696, 0.963, 0.965, 0.525, 0.876, 0.375, 0.767, 0.414, 0.942, 0.98, 0.974, 0.968, 0.974, 0.993, 0.607, 0.957, 0.832, 0.399, 0.771, 0.839, 0.831, 0.815, 0.46, 0.964, 0.914, 0.456, 0.895, 0.739, 0.898, 0.959, 0.673, 0.876, 0.988, 0.79, 0.821, 0.598, 0.922, 0.84, 0.915, 0.995, 0.622, 0.998, 0.936, 0.963, 0.969, 0.637, 0.882, 0.901, 0.477, 0.708, 0.919, 0.783, 0.952, 0.181, 0.782, 0.514, 0.946, 0.929, 0.659, 0.852, 0.969, 0.453, 0.738, 0.927, 0.959, 0.705, 0.495, 0.927, 0.824, 0.566, 0.965, 0.513]
global origin = 1
global destination = 40