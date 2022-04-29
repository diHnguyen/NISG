global arcs = [1 11; 1 15; 2 7; 2 23; 2 27; 2 33; 2 37; 3 16; 4 7; 4 10; 4 13; 4 36; 5 11; 5 15; 5 26; 5 36; 6 17; 6 30; 7 14; 7 20; 7 35; 7 40; 8 14; 8 15; 9 3; 9 21; 10 3; 10 8; 10 17; 11 5; 11 15; 11 23; 11 31; 11 32; 11 33; 12 22; 13 4; 13 32; 14 11; 14 19; 14 22; 14 25; 14 33; 14 39; 15 23; 15 24; 15 26; 16 33; 16 34; 17 26; 17 37; 18 5; 18 29; 19 3; 19 8; 19 17; 19 25; 19 40; 20 18; 20 29; 20 32; 20 35; 21 5; 21 8; 21 17; 21 30; 21 33; 22 27; 22 28; 22 31; 23 4; 23 7; 23 10; 23 24; 23 31; 24 2; 24 19; 24 33; 24 38; 25 12; 26 2; 26 6; 26 8; 26 10; 27 4; 27 5; 27 9; 27 14; 27 28; 28 19; 28 29; 28 35; 29 7; 29 9; 29 13; 29 22; 29 33; 30 17; 30 20; 30 26; 30 36; 30 40; 31 10; 31 15; 31 21; 32 13; 32 14; 33 10; 33 11; 33 30; 34 21; 34 35; 34 38; 34 39; 35 5; 35 13; 35 17; 35 24; 36 2; 36 9; 36 20; 36 37; 37 9; 38 11; 38 14; 38 18; 38 24; 39 27; 39 31; 39 38]
global d_x = [8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 10.0, 8.0, 10.0, 10.0, 1.0, 2.0, 6.0, 1.0, 9.0, 8.0, 4.0, 3.0, 9.0, 3.0, 3.0, 2.0, 10.0, 6.0, 2.0, 9.0, 3.0, 7.0, 2.0, 7.0, 8.0, 10.0, 4.0, 8.0, 7.0, 6.0, 8.0, 6.0, 8.0, 8.0, 7.0, 2.0, 3.0, 5.0, 1.0, 10.0, 5.0, 8.0, 10.0, 9.0, 5.0, 6.0, 8.0, 2.0, 6.0, 10.0, 2.0, 4.0, 5.0, 6.0, 7.0, 2.0, 3.0, 1.0, 3.0, 2.0, 10.0, 8.0, 2.0, 6.0, 4.0, 9.0, 1.0, 9.0, 3.0, 6.0, 2.0, 10.0, 2.0, 10.0, 7.0, 8.0, 4.0, 8.0, 2.0, 1.0, 3.0, 5.0, 4.0, 4.0, 1.0, 5.0, 4.0, 4.0, 6.0, 10.0, 5.0, 1.0, 5.0, 5.0, 5.0, 10.0, 7.0, 1.0, 5.0, 6.0, 1.0, 5.0, 8.0, 10.0, 3.0, 8.0, 10.0, 2.0, 9.0, 6.0, 3.0, 5.0, 4.0, 3.0, 6.0, 8.0, 5.0, 3.0, 6.0, 9.0, 1.0, 3.0, 9.0, 9.0]
global b_x = 5
global d_y = [1.0, 3.0, 2.0, 8.0, 6.0, 9.0, 8.0, 8.0, 3.0, 10.0, 4.0, 1.0, 8.0, 3.0, 9.0, 9.0, 6.0, 8.0, 7.0, 5.0, 8.0, 7.0, 6.0, 1.0, 7.0, 8.0, 3.0, 6.0, 8.0, 8.0, 5.0, 2.0, 3.0, 8.0, 6.0, 9.0, 1.0, 2.0, 3.0, 1.0, 5.0, 1.0, 8.0, 4.0, 7.0, 10.0, 3.0, 2.0, 4.0, 3.0, 1.0, 1.0, 7.0, 9.0, 2.0, 5.0, 6.0, 2.0, 5.0, 7.0, 1.0, 10.0, 8.0, 1.0, 1.0, 7.0, 10.0, 1.0, 10.0, 1.0, 7.0, 7.0, 8.0, 1.0, 2.0, 9.0, 5.0, 10.0, 2.0, 8.0, 4.0, 2.0, 6.0, 7.0, 6.0, 1.0, 2.0, 10.0, 10.0, 3.0, 5.0, 9.0, 9.0, 10.0, 3.0, 10.0, 10.0, 9.0, 5.0, 9.0, 2.0, 5.0, 3.0, 3.0, 9.0, 7.0, 7.0, 5.0, 5.0, 5.0, 1.0, 3.0, 4.0, 3.0, 5.0, 9.0, 6.0, 9.0, 2.0, 7.0, 3.0, 5.0, 10.0, 4.0, 3.0, 8.0, 8.0, 7.0, 3.0, 10.0]
global b_y = 10
global p = [0.006, 0.653, 0.91, 0.497, 0.966, 0.381, 0.174, 0.207, 0.629, 0.562, 0.282, 0.323, 0.049, 0.361, 0.562, 0.931, 0.775, 0.493, 0.006, 0.993, 0.84, 0.694, 0.167, 0.064, 0.351, 0.148, 0.284, 0.304, 0.843, 0.661, 0.493, 0.8, 0.865, 0.039, 0.191, 0.054, 0.106, 0.36, 0.502, 0.49, 0.249, 0.937, 0.81, 0.774, 0.36, 0.5, 0.863, 0.853, 0.41, 0.592, 0.102, 0.487, 0.548, 0.355, 0.145, 0.467, 0.085, 0.88, 0.064, 0.17, 0.423, 0.709, 0.262, 0.874, 0.401, 0.929, 0.424, 0.379, 0.171, 0.154, 0.15, 0.823, 0.031, 0.417, 0.121, 0.5, 0.0, 0.974, 0.943, 0.621, 0.172, 0.496, 0.43, 0.638, 0.945, 0.393, 0.818, 0.778, 0.235, 0.43, 0.347, 0.319, 0.266, 0.91, 0.92, 0.089, 0.587, 0.031, 0.654, 0.088, 0.164, 0.068, 0.271, 0.092, 0.197, 0.337, 0.902, 0.481, 0.783, 0.547, 0.583, 0.7, 0.132, 0.134, 0.542, 0.086, 0.607, 0.647, 0.137, 0.916, 0.006, 0.252, 0.62, 0.485, 0.205, 0.259, 0.165, 0.932, 0.615, 0.473]
global q = [0.369, 0.813, 0.998, 0.811, 0.973, 0.981, 0.339, 0.656, 0.759, 0.643, 0.659, 0.341, 0.137, 0.483, 0.759, 0.932, 0.919, 0.762, 0.693, 0.995, 0.862, 0.959, 0.494, 0.691, 0.411, 0.734, 0.887, 0.323, 0.967, 0.742, 0.816, 0.842, 0.914, 0.192, 0.391, 0.269, 0.857, 0.923, 0.866, 0.621, 0.642, 0.97, 0.835, 0.995, 0.652, 0.677, 0.982, 0.923, 0.53, 0.794, 0.993, 0.557, 0.686, 0.611, 0.329, 0.879, 0.328, 0.982, 0.533, 0.366, 0.427, 0.78, 0.666, 0.911, 0.472, 0.99, 0.497, 0.896, 0.477, 0.538, 0.168, 0.955, 0.401, 0.624, 0.776, 0.505, 0.922, 0.993, 0.981, 0.818, 0.419, 0.958, 0.902, 0.819, 0.966, 0.575, 0.828, 0.914, 0.25, 0.45, 0.657, 0.942, 0.396, 0.915, 0.982, 0.413, 0.871, 0.474, 0.946, 0.513, 0.775, 0.749, 0.557, 0.436, 0.421, 0.528, 0.908, 0.88, 0.963, 0.744, 0.961, 0.853, 0.625, 0.778, 0.567, 0.783, 0.648, 0.692, 0.909, 0.963, 0.138, 0.631, 0.872, 0.59, 0.422, 0.685, 0.903, 0.951, 0.853, 0.991]
global origin = 1
global destination = 40
