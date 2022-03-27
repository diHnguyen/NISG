global arcs = [1 5; 1 16; 1 33; 1 39; 1 40; 1 41; 1 50; 2 31; 2 36; 2 39; 2 45; 2 50; 3 18; 3 41; 3 45; 3 49; 4 13; 4 14; 4 22; 4 24; 4 32; 4 42; 4 43; 5 2; 5 17; 5 22; 5 32; 5 43; 5 50; 6 14; 6 21; 6 38; 6 44; 7 42; 8 5; 8 24; 8 44; 8 48; 9 8; 9 13; 9 22; 9 23; 9 28; 9 30; 9 33; 9 37; 9 44; 10 4; 10 19; 10 29; 10 43; 11 3; 11 16; 11 25; 11 29; 11 32; 11 35; 11 37; 11 48; 12 40; 12 46; 13 16; 13 19; 13 28; 13 30; 13 48; 14 34; 15 13; 15 25; 15 29; 15 36; 16 4; 16 7; 16 22; 17 33; 17 42; 18 4; 18 7; 18 27; 18 32; 18 44; 19 3; 19 6; 19 10; 19 11; 19 26; 19 37; 19 48; 20 3; 20 9; 20 13; 20 27; 20 31; 21 2; 21 6; 21 7; 21 28; 21 34; 21 40; 22 9; 23 12; 23 13; 23 29; 23 36; 23 41; 24 13; 24 20; 24 34; 24 42; 25 12; 25 34; 25 50; 26 8; 26 21; 26 24; 26 44; 26 47; 27 6; 27 7; 27 42; 28 4; 28 5; 28 6; 28 21; 28 31; 28 33; 28 37; 28 39; 28 41; 29 10; 29 12; 29 19; 29 36; 29 42; 30 9; 30 26; 30 40; 31 5; 31 20; 31 32; 31 38; 31 40; 32 12; 32 35; 32 46; 33 4; 33 11; 33 14; 33 20; 34 2; 34 12; 34 14; 35 3; 35 8; 35 24; 35 36; 36 4; 37 15; 37 20; 37 25; 38 8; 38 11; 38 25; 38 33; 38 40; 38 45; 39 22; 39 24; 39 27; 39 42; 39 44; 40 4; 40 5; 40 7; 40 9; 40 17; 40 32; 40 34; 40 44; 41 6; 41 12; 41 13; 41 21; 41 27; 41 34; 42 26; 42 32; 42 41; 42 49; 43 4; 43 37; 44 6; 44 26; 44 40; 45 7; 45 16; 45 20; 45 29; 45 31; 45 49; 46 7; 46 8; 46 15; 46 16; 46 22; 46 28; 46 36; 46 39; 47 3; 47 16; 47 24; 47 26; 47 31; 47 37; 47 40; 48 2; 48 4; 48 8; 48 19; 48 35; 48 45; 48 50; 49 12; 49 20; 49 24; 49 36; 49 47]
global d_x = [10.0, 6.0, 6.0, 7.0, 6.0, 5.0, 2.0, 10.0, 5.0, 10.0, 10.0, 6.0, 7.0, 7.0, 2.0, 8.0, 3.0, 1.0, 3.0, 4.0, 4.0, 7.0, 7.0, 9.0, 7.0, 5.0, 2.0, 10.0, 10.0, 2.0, 5.0, 1.0, 2.0, 10.0, 8.0, 10.0, 1.0, 7.0, 3.0, 2.0, 5.0, 5.0, 2.0, 8.0, 9.0, 2.0, 2.0, 2.0, 4.0, 8.0, 5.0, 7.0, 5.0, 1.0, 5.0, 10.0, 9.0, 8.0, 3.0, 5.0, 1.0, 10.0, 3.0, 7.0, 9.0, 7.0, 7.0, 6.0, 1.0, 9.0, 10.0, 6.0, 9.0, 5.0, 3.0, 7.0, 1.0, 7.0, 1.0, 2.0, 9.0, 7.0, 10.0, 1.0, 9.0, 8.0, 8.0, 8.0, 10.0, 3.0, 9.0, 6.0, 9.0, 3.0, 1.0, 9.0, 5.0, 2.0, 5.0, 5.0, 10.0, 4.0, 1.0, 9.0, 3.0, 5.0, 9.0, 6.0, 7.0, 10.0, 3.0, 1.0, 3.0, 2.0, 8.0, 6.0, 2.0, 10.0, 7.0, 9.0, 6.0, 7.0, 4.0, 4.0, 2.0, 6.0, 10.0, 6.0, 5.0, 8.0, 1.0, 8.0, 9.0, 7.0, 2.0, 1.0, 6.0, 10.0, 10.0, 10.0, 8.0, 10.0, 3.0, 2.0, 6.0, 9.0, 5.0, 1.0, 10.0, 1.0, 3.0, 1.0, 5.0, 3.0, 10.0, 7.0, 3.0, 6.0, 2.0, 2.0, 7.0, 1.0, 3.0, 6.0, 5.0, 2.0, 6.0, 6.0, 1.0, 8.0, 10.0, 10.0, 9.0, 10.0, 10.0, 8.0, 1.0, 8.0, 8.0, 3.0, 8.0, 10.0, 6.0, 3.0, 10.0, 6.0, 9.0, 5.0, 6.0, 2.0, 4.0, 7.0, 6.0, 5.0, 3.0, 1.0, 10.0, 2.0, 3.0, 2.0, 1.0, 6.0, 6.0, 8.0, 7.0, 6.0, 10.0, 3.0, 2.0, 9.0, 7.0, 3.0, 5.0, 2.0, 7.0, 7.0, 1.0, 5.0, 7.0, 6.0, 3.0, 8.0, 5.0, 5.0, 1.0, 4.0, 8.0]
global b_x = 5
global d_y = [5.0, 9.0, 2.0, 7.0, 7.0, 8.0, 10.0, 5.0, 2.0, 5.0, 8.0, 6.0, 4.0, 8.0, 7.0, 7.0, 7.0, 7.0, 9.0, 7.0, 9.0, 7.0, 10.0, 10.0, 10.0, 7.0, 10.0, 6.0, 1.0, 9.0, 10.0, 10.0, 4.0, 6.0, 8.0, 1.0, 10.0, 10.0, 8.0, 9.0, 2.0, 6.0, 8.0, 1.0, 5.0, 1.0, 6.0, 2.0, 9.0, 2.0, 10.0, 10.0, 7.0, 8.0, 7.0, 10.0, 3.0, 10.0, 9.0, 10.0, 6.0, 10.0, 5.0, 4.0, 4.0, 10.0, 9.0, 2.0, 1.0, 8.0, 6.0, 2.0, 9.0, 2.0, 8.0, 7.0, 6.0, 7.0, 6.0, 4.0, 5.0, 1.0, 5.0, 3.0, 10.0, 2.0, 9.0, 3.0, 4.0, 6.0, 8.0, 6.0, 2.0, 10.0, 7.0, 4.0, 9.0, 3.0, 5.0, 10.0, 4.0, 5.0, 5.0, 9.0, 4.0, 1.0, 6.0, 4.0, 10.0, 5.0, 6.0, 5.0, 8.0, 8.0, 7.0, 10.0, 2.0, 5.0, 6.0, 6.0, 9.0, 3.0, 4.0, 1.0, 10.0, 8.0, 3.0, 1.0, 1.0, 3.0, 5.0, 5.0, 10.0, 7.0, 5.0, 8.0, 9.0, 3.0, 2.0, 5.0, 10.0, 7.0, 8.0, 1.0, 9.0, 8.0, 10.0, 7.0, 10.0, 8.0, 4.0, 3.0, 9.0, 7.0, 3.0, 4.0, 1.0, 4.0, 7.0, 4.0, 3.0, 2.0, 7.0, 10.0, 2.0, 10.0, 6.0, 6.0, 2.0, 3.0, 9.0, 3.0, 9.0, 1.0, 8.0, 4.0, 1.0, 3.0, 5.0, 10.0, 4.0, 8.0, 7.0, 2.0, 9.0, 1.0, 6.0, 2.0, 10.0, 9.0, 2.0, 1.0, 5.0, 10.0, 2.0, 9.0, 1.0, 9.0, 2.0, 5.0, 7.0, 1.0, 10.0, 8.0, 2.0, 7.0, 3.0, 2.0, 5.0, 5.0, 1.0, 1.0, 10.0, 9.0, 1.0, 5.0, 8.0, 9.0, 4.0, 8.0, 6.0, 1.0, 7.0, 2.0, 7.0, 6.0, 4.0]
global b_y = 10
global p = [0.8, 0.54, 0.444, 0.493, 0.888, 0.152, 0.69, 0.307, 0.5, 0.334, 0.654, 0.842, 0.836, 0.519, 0.805, 0.069, 0.99, 0.527, 0.876, 0.742, 0.456, 0.959, 0.359, 0.69, 0.838, 0.201, 0.05, 0.3, 0.627, 0.792, 0.854, 0.877, 0.548, 0.878, 0.052, 0.936, 0.306, 0.826, 0.266, 0.382, 0.076, 0.366, 0.408, 0.89, 0.836, 0.479, 0.413, 0.006, 0.985, 0.177, 0.465, 0.984, 0.316, 0.695, 0.483, 0.339, 0.643, 0.502, 0.462, 0.757, 0.897, 0.853, 0.143, 0.993, 0.596, 0.739, 0.792, 0.361, 0.361, 0.736, 0.927, 0.894, 0.731, 0.359, 0.563, 0.745, 0.431, 0.936, 0.944, 0.359, 0.671, 0.328, 0.654, 0.395, 0.032, 0.416, 0.117, 0.316, 0.379, 0.368, 0.718, 0.946, 0.953, 0.665, 0.745, 0.795, 0.946, 0.245, 0.171, 0.095, 0.989, 0.953, 0.735, 0.765, 0.032, 0.599, 0.782, 0.856, 0.556, 0.623, 0.987, 0.894, 0.537, 0.612, 0.725, 0.923, 0.654, 0.168, 0.716, 0.086, 0.681, 0.52, 0.078, 0.133, 0.129, 0.748, 0.803, 0.682, 0.008, 0.557, 0.46, 0.485, 0.188, 0.988, 0.69, 0.288, 0.893, 0.589, 0.237, 0.006, 0.786, 0.932, 0.839, 0.738, 0.271, 0.886, 0.591, 0.758, 0.795, 0.229, 0.027, 0.753, 0.491, 0.981, 0.805, 0.271, 0.615, 0.143, 0.947, 0.606, 0.83, 0.61, 0.059, 0.955, 0.132, 0.858, 0.918, 0.9, 0.95, 0.437, 0.385, 0.862, 0.909, 0.264, 0.756, 0.966, 0.715, 0.716, 0.698, 0.022, 0.224, 0.829, 0.172, 0.66, 0.35, 0.965, 0.196, 0.683, 0.676, 0.377, 0.265, 0.21, 0.264, 0.204, 0.307, 0.338, 0.983, 0.616, 0.792, 0.505, 0.847, 0.783, 0.65, 0.515, 0.376, 0.714, 0.52, 0.578, 0.246, 0.117, 0.045, 0.358, 0.188, 0.325, 0.322, 0.315, 0.041, 0.764, 0.694, 0.163, 0.134, 0.564, 0.104, 0.762, 0.711, 0.997, 0.959]
global q = [0.984, 0.804, 0.698, 0.593, 0.978, 0.56, 0.991, 0.854, 0.822, 0.929, 0.685, 0.919, 0.964, 0.612, 0.989, 0.705, 0.999, 0.56, 0.926, 0.985, 0.505, 0.99, 0.487, 0.897, 0.908, 0.948, 0.125, 0.66, 0.973, 0.906, 0.877, 0.998, 0.652, 0.944, 0.321, 0.988, 0.395, 0.964, 0.549, 0.647, 0.169, 0.994, 0.583, 0.918, 0.991, 0.517, 0.497, 0.093, 0.995, 0.844, 0.781, 0.992, 0.666, 0.918, 0.896, 0.998, 0.904, 0.565, 0.949, 0.927, 0.917, 0.999, 0.197, 0.994, 0.906, 0.742, 0.887, 0.756, 0.729, 0.821, 0.962, 0.92, 0.906, 0.583, 0.849, 0.772, 0.478, 0.94, 0.951, 0.69, 0.79, 0.456, 0.867, 0.661, 0.976, 0.732, 0.189, 0.713, 0.978, 0.824, 0.96, 0.955, 0.977, 0.836, 0.994, 0.799, 0.997, 0.82, 0.187, 0.888, 0.997, 0.984, 0.936, 0.942, 0.88, 0.871, 0.94, 0.922, 0.911, 0.837, 0.995, 0.987, 0.812, 0.659, 0.737, 0.97, 0.938, 0.656, 0.927, 0.587, 0.895, 0.902, 0.829, 0.28, 0.296, 0.906, 0.878, 0.788, 0.221, 0.623, 0.818, 0.952, 0.531, 0.993, 0.768, 0.556, 0.948, 0.93, 0.513, 0.335, 0.808, 0.962, 0.951, 0.948, 0.707, 0.891, 0.939, 0.803, 0.879, 0.625, 0.054, 0.869, 0.84, 0.997, 0.978, 0.695, 0.748, 0.299, 0.962, 0.642, 0.858, 0.825, 0.568, 0.975, 0.863, 0.981, 0.925, 0.935, 0.982, 0.883, 0.533, 0.898, 0.932, 0.664, 0.789, 0.987, 0.881, 0.93, 0.763, 0.532, 0.56, 0.961, 0.804, 0.848, 0.492, 0.966, 0.922, 0.962, 0.983, 0.62, 0.435, 0.853, 0.384, 0.22, 0.89, 0.669, 0.988, 0.905, 0.845, 0.653, 0.978, 0.985, 0.897, 0.93, 0.924, 0.909, 0.612, 0.802, 0.774, 0.28, 0.296, 0.691, 0.993, 0.59, 0.725, 0.878, 0.343, 0.902, 0.771, 0.176, 0.211, 0.995, 0.356, 0.94, 0.879, 0.998, 0.974]
global origin = 1
global destination = 50