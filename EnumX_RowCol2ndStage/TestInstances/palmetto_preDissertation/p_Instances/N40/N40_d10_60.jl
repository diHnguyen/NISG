global arcs = [1 25; 1 37; 2 7; 2 13; 2 21; 2 28; 2 36; 2 39; 3 5; 3 6; 3 23; 3 33; 4 16; 4 22; 5 8; 5 17; 5 20; 5 32; 6 3; 6 5; 6 11; 6 22; 6 30; 6 39; 7 18; 7 26; 7 30; 7 40; 8 4; 8 6; 8 12; 8 17; 8 25; 8 30; 8 39; 9 12; 9 24; 9 30; 9 39; 10 6; 10 23; 10 36; 11 15; 11 20; 11 23; 11 26; 11 27; 11 31; 11 33; 11 37; 11 40; 12 30; 13 8; 13 17; 13 21; 13 34; 13 36; 13 39; 14 12; 14 13; 14 20; 14 35; 15 3; 15 7; 15 8; 15 18; 15 19; 15 29; 15 34; 15 37; 16 4; 16 7; 16 35; 17 10; 17 12; 17 38; 18 19; 18 36; 19 3; 19 4; 19 5; 19 16; 19 37; 20 3; 20 4; 21 2; 21 37; 22 6; 22 17; 22 18; 22 23; 23 3; 23 22; 23 27; 24 4; 24 19; 25 15; 25 28; 25 40; 26 8; 26 10; 26 23; 27 6; 27 7; 27 21; 28 3; 28 9; 28 19; 28 36; 29 4; 29 39; 30 5; 30 10; 30 12; 30 26; 31 22; 32 23; 32 26; 32 28; 32 33; 32 38; 33 21; 33 31; 33 32; 34 3; 34 6; 34 10; 34 14; 34 23; 34 39; 35 3; 35 34; 36 29; 37 9; 37 19; 37 20; 37 21; 37 30; 37 31; 38 27; 38 35; 39 6]
global d_x = [5.0, 2.0, 2.0, 4.0, 4.0, 4.0, 5.0, 6.0, 8.0, 1.0, 8.0, 10.0, 4.0, 1.0, 3.0, 10.0, 8.0, 1.0, 7.0, 9.0, 9.0, 8.0, 10.0, 1.0, 2.0, 1.0, 8.0, 6.0, 4.0, 6.0, 4.0, 6.0, 3.0, 10.0, 7.0, 1.0, 8.0, 7.0, 10.0, 9.0, 3.0, 7.0, 10.0, 6.0, 5.0, 2.0, 5.0, 6.0, 3.0, 2.0, 4.0, 5.0, 10.0, 7.0, 5.0, 1.0, 9.0, 5.0, 4.0, 9.0, 4.0, 5.0, 5.0, 2.0, 10.0, 1.0, 8.0, 6.0, 2.0, 10.0, 2.0, 6.0, 6.0, 10.0, 5.0, 10.0, 9.0, 5.0, 5.0, 3.0, 1.0, 6.0, 6.0, 9.0, 10.0, 3.0, 6.0, 3.0, 9.0, 5.0, 6.0, 7.0, 5.0, 1.0, 1.0, 4.0, 5.0, 3.0, 5.0, 9.0, 1.0, 3.0, 4.0, 5.0, 5.0, 1.0, 5.0, 6.0, 8.0, 4.0, 4.0, 9.0, 3.0, 3.0, 3.0, 2.0, 2.0, 4.0, 2.0, 9.0, 2.0, 1.0, 4.0, 6.0, 9.0, 5.0, 2.0, 7.0, 5.0, 2.0, 5.0, 4.0, 6.0, 6.0, 8.0, 9.0, 3.0, 7.0, 7.0, 4.0, 4.0, 2.0]
global b_x = 5
global d_y = [9.0, 1.0, 7.0, 4.0, 2.0, 3.0, 4.0, 4.0, 3.0, 4.0, 2.0, 5.0, 4.0, 3.0, 5.0, 6.0, 7.0, 2.0, 7.0, 4.0, 1.0, 4.0, 4.0, 1.0, 8.0, 6.0, 4.0, 4.0, 1.0, 3.0, 5.0, 9.0, 6.0, 2.0, 4.0, 10.0, 1.0, 5.0, 10.0, 6.0, 7.0, 1.0, 6.0, 4.0, 2.0, 10.0, 1.0, 7.0, 9.0, 3.0, 6.0, 6.0, 8.0, 3.0, 3.0, 9.0, 1.0, 6.0, 6.0, 9.0, 6.0, 5.0, 4.0, 3.0, 5.0, 7.0, 6.0, 8.0, 5.0, 9.0, 4.0, 10.0, 7.0, 5.0, 6.0, 4.0, 9.0, 4.0, 4.0, 4.0, 8.0, 1.0, 2.0, 5.0, 4.0, 6.0, 8.0, 2.0, 10.0, 5.0, 7.0, 4.0, 3.0, 4.0, 4.0, 4.0, 3.0, 1.0, 4.0, 6.0, 5.0, 7.0, 7.0, 10.0, 7.0, 6.0, 6.0, 2.0, 1.0, 10.0, 8.0, 4.0, 10.0, 1.0, 3.0, 6.0, 6.0, 9.0, 10.0, 10.0, 4.0, 5.0, 9.0, 9.0, 5.0, 2.0, 1.0, 3.0, 2.0, 7.0, 1.0, 8.0, 2.0, 1.0, 9.0, 10.0, 8.0, 2.0, 7.0, 3.0, 7.0, 5.0]
global b_y = 10
global p = [0.654, 0.536, 0.935, 0.066, 0.024, 0.572, 0.485, 0.972, 0.779, 0.632, 0.831, 0.08, 0.639, 0.28, 0.586, 0.938, 0.558, 0.221, 0.403, 0.218, 0.24, 0.547, 0.952, 0.434, 0.992, 0.295, 0.967, 0.379, 0.607, 0.392, 0.708, 0.354, 0.028, 0.812, 0.836, 0.782, 0.006, 0.149, 0.027, 0.488, 0.488, 0.07, 0.104, 0.204, 0.901, 0.932, 0.908, 0.586, 0.15, 0.784, 0.55, 0.355, 0.694, 0.331, 0.485, 0.476, 0.433, 0.312, 0.638, 0.72, 0.577, 0.172, 0.501, 0.271, 0.075, 0.37, 0.958, 0.798, 0.683, 0.455, 0.591, 0.06, 0.629, 0.527, 0.627, 0.172, 0.23, 0.264, 0.158, 0.032, 0.549, 0.07, 0.478, 0.7, 0.891, 0.219, 0.128, 0.652, 0.835, 0.933, 0.802, 0.665, 0.074, 0.477, 0.297, 0.708, 0.939, 0.558, 0.143, 0.306, 0.99, 0.775, 0.083, 0.299, 0.434, 0.138, 0.052, 0.582, 0.832, 0.773, 0.842, 0.447, 0.708, 0.139, 0.986, 0.797, 0.457, 0.748, 0.299, 0.829, 0.612, 0.877, 0.479, 0.486, 0.943, 0.161, 0.333, 0.883, 0.33, 0.061, 0.225, 0.67, 0.401, 0.072, 0.896, 0.764, 0.751, 0.453, 0.656, 0.12, 0.637, 0.73]
global q = [0.962, 0.699, 0.987, 0.123, 0.491, 0.743, 0.894, 0.997, 0.944, 0.902, 0.855, 0.312, 0.982, 0.964, 0.806, 0.991, 0.755, 0.785, 0.484, 0.979, 0.55, 0.783, 0.971, 0.56, 0.993, 0.795, 0.996, 0.939, 0.785, 0.921, 0.756, 0.383, 0.165, 0.836, 0.872, 0.919, 0.347, 0.56, 0.165, 0.707, 0.723, 0.308, 0.227, 0.406, 0.926, 0.976, 0.999, 0.83, 0.836, 0.953, 0.781, 0.857, 0.853, 0.469, 0.903, 0.733, 0.556, 0.638, 0.757, 0.861, 0.965, 0.615, 0.748, 0.983, 0.973, 0.861, 0.981, 0.99, 0.705, 0.738, 0.656, 0.591, 0.716, 0.583, 0.694, 0.624, 0.764, 0.652, 0.393, 0.586, 0.656, 0.716, 0.948, 0.804, 0.987, 0.787, 0.998, 0.834, 0.936, 0.998, 0.992, 0.898, 0.611, 0.58, 0.926, 0.927, 0.967, 0.713, 0.183, 0.74, 0.998, 0.928, 0.489, 0.33, 0.819, 0.633, 0.314, 0.714, 0.964, 0.977, 0.888, 0.969, 0.999, 0.342, 0.995, 0.906, 0.762, 0.769, 0.531, 0.95, 0.774, 0.877, 0.984, 0.666, 0.967, 0.324, 0.952, 0.966, 0.821, 0.09, 0.274, 0.845, 0.747, 0.343, 0.902, 0.947, 0.794, 0.901, 0.667, 0.195, 0.934, 0.871]
global origin = 1
global destination = 40