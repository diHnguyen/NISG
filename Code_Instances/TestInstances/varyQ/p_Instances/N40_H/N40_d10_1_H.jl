global arcs = [1 10; 1 32; 2 9; 2 29; 2 39; 3 7; 3 13; 3 16; 4 3; 4 21; 4 26; 5 19; 5 20; 5 26; 5 30; 6 7; 6 8; 6 14; 6 27; 6 30; 6 40; 7 12; 7 13; 7 23; 7 39; 8 3; 8 18; 9 38; 10 29; 10 30; 11 8; 11 17; 11 23; 11 25; 11 27; 11 32; 11 35; 12 2; 12 13; 12 16; 12 31; 12 32; 13 14; 13 16; 13 31; 13 34; 14 16; 14 27; 14 38; 15 8; 15 14; 15 25; 15 28; 16 9; 16 12; 16 21; 16 29; 16 32; 16 33; 17 4; 17 12; 17 27; 17 32; 17 38; 18 22; 18 25; 18 28; 18 35; 18 39; 19 3; 19 7; 19 12; 19 33; 20 13; 20 18; 20 24; 20 28; 20 33; 20 36; 20 39; 21 3; 21 12; 21 13; 21 16; 22 7; 22 11; 22 12; 22 13; 22 26; 22 33; 22 38; 22 39; 23 15; 23 22; 23 30; 23 39; 24 6; 24 10; 24 20; 24 22; 24 23; 24 25; 24 30; 24 35; 25 16; 25 24; 26 4; 26 24; 27 6; 27 15; 27 31; 27 35; 27 40; 28 14; 28 15; 28 24; 29 8; 29 16; 29 18; 29 19; 29 35; 30 5; 30 7; 30 22; 30 31; 30 37; 31 14; 31 26; 31 40; 32 5; 32 7; 32 12; 32 17; 32 21; 32 30; 33 5; 33 7; 33 10; 33 14; 33 29; 34 5; 34 22; 34 28; 35 27; 35 30; 35 33; 36 4; 36 12; 36 13; 36 17; 37 9; 37 10; 37 25; 37 29; 38 21; 38 27; 39 2; 39 7; 39 15; 39 22; 39 28; 39 31; 39 32]
global d_x = [5.0, 2.0, 2.0, 5.0, 4.0, 4.0, 2.0, 6.0, 3.0, 5.0, 6.0, 7.0, 7.0, 9.0, 1.0, 9.0, 9.0, 6.0, 4.0, 1.0, 4.0, 6.0, 5.0, 10.0, 5.0, 3.0, 2.0, 7.0, 7.0, 2.0, 10.0, 8.0, 2.0, 8.0, 8.0, 1.0, 6.0, 10.0, 3.0, 5.0, 3.0, 6.0, 1.0, 7.0, 4.0, 8.0, 2.0, 2.0, 6.0, 7.0, 2.0, 1.0, 8.0, 6.0, 1.0, 5.0, 5.0, 9.0, 2.0, 6.0, 1.0, 2.0, 6.0, 7.0, 4.0, 4.0, 4.0, 8.0, 7.0, 3.0, 7.0, 9.0, 3.0, 8.0, 8.0, 2.0, 9.0, 4.0, 5.0, 8.0, 3.0, 5.0, 7.0, 7.0, 9.0, 8.0, 1.0, 9.0, 1.0, 9.0, 1.0, 7.0, 6.0, 8.0, 6.0, 7.0, 9.0, 9.0, 2.0, 1.0, 3.0, 1.0, 1.0, 9.0, 7.0, 2.0, 5.0, 3.0, 6.0, 1.0, 10.0, 10.0, 1.0, 7.0, 1.0, 4.0, 4.0, 6.0, 4.0, 5.0, 3.0, 7.0, 8.0, 5.0, 9.0, 6.0, 6.0, 7.0, 3.0, 2.0, 3.0, 7.0, 6.0, 7.0, 8.0, 4.0, 8.0, 3.0, 9.0, 1.0, 8.0, 6.0, 1.0, 3.0, 1.0, 1.0, 2.0, 9.0, 7.0, 9.0, 1.0, 9.0, 4.0, 9.0, 7.0, 8.0, 1.0, 9.0, 10.0, 4.0, 5.0, 10.0, 4.0]
global b_x = 5
global d_y = [3.0, 1.0, 9.0, 7.0, 10.0, 1.0, 7.0, 5.0, 3.0, 10.0, 1.0, 2.0, 6.0, 1.0, 2.0, 9.0, 5.0, 5.0, 6.0, 9.0, 5.0, 4.0, 8.0, 5.0, 4.0, 9.0, 9.0, 1.0, 8.0, 6.0, 9.0, 6.0, 7.0, 5.0, 4.0, 6.0, 2.0, 10.0, 10.0, 1.0, 5.0, 1.0, 4.0, 6.0, 5.0, 5.0, 7.0, 1.0, 6.0, 1.0, 6.0, 9.0, 1.0, 7.0, 9.0, 4.0, 7.0, 9.0, 4.0, 6.0, 2.0, 10.0, 3.0, 8.0, 1.0, 3.0, 2.0, 5.0, 7.0, 5.0, 8.0, 7.0, 1.0, 10.0, 7.0, 8.0, 10.0, 10.0, 9.0, 5.0, 9.0, 9.0, 6.0, 2.0, 9.0, 1.0, 4.0, 8.0, 9.0, 1.0, 8.0, 7.0, 4.0, 7.0, 9.0, 2.0, 2.0, 6.0, 2.0, 3.0, 2.0, 4.0, 9.0, 7.0, 1.0, 9.0, 5.0, 3.0, 2.0, 5.0, 6.0, 1.0, 8.0, 3.0, 4.0, 4.0, 5.0, 4.0, 1.0, 10.0, 8.0, 4.0, 3.0, 1.0, 3.0, 7.0, 9.0, 10.0, 1.0, 10.0, 9.0, 5.0, 1.0, 9.0, 3.0, 1.0, 7.0, 10.0, 1.0, 1.0, 2.0, 6.0, 8.0, 10.0, 5.0, 1.0, 4.0, 6.0, 2.0, 3.0, 1.0, 9.0, 2.0, 3.0, 8.0, 9.0, 6.0, 4.0, 8.0, 7.0, 9.0, 10.0, 4.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.853, 0.927, 0.962, 0.864, 0.985, 0.967, 0.894, 0.905, 0.838, 0.816, 0.811, 0.816, 0.924, 0.964, 0.853, 0.917, 0.985, 0.974, 0.948, 0.836, 0.825, 0.956, 0.931, 0.958, 0.862, 0.89, 0.835, 0.992, 0.86, 0.957, 0.877, 0.963, 0.848, 0.943, 0.98, 0.88, 0.95, 0.949, 0.807, 0.863, 0.943, 0.957, 0.885, 0.854, 0.973, 0.913, 0.975, 0.998, 0.863, 0.97, 0.938, 0.903, 0.834, 0.877, 0.907, 0.98, 0.956, 0.901, 0.979, 0.988, 0.957, 0.98, 0.957, 0.99, 0.925, 0.906, 0.859, 0.848, 0.802, 0.963, 0.933, 0.963, 0.957, 0.838, 0.853, 0.841, 0.879, 0.871, 0.823, 0.877, 0.837, 0.925, 0.873, 0.861, 0.828, 0.815, 0.813, 0.935, 0.998, 0.884, 0.888, 0.954, 0.847, 0.963, 0.967, 0.966, 0.989, 0.9, 0.857, 0.808, 0.91, 0.996, 0.906, 0.913, 0.873, 0.98, 0.81, 0.854, 0.864, 0.935, 0.977, 0.835, 0.886, 0.805, 0.828, 0.952, 0.936, 0.966, 0.92, 0.875, 0.88, 0.926, 0.821, 0.817, 0.808, 0.874, 0.81, 0.846, 0.917, 0.884, 0.851, 0.969, 0.83, 0.891, 0.824, 0.803, 0.838, 0.846, 0.944, 0.999, 0.939, 0.835, 0.981, 0.818, 0.85, 0.981, 0.89, 0.911, 0.857, 0.926, 0.973, 0.858, 0.863, 0.993, 0.928, 0.986, 0.834, 0.879, 0.997, 0.936, 0.876, 0.901, 0.909]
global origin = 1
global destination = 40