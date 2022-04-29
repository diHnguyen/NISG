global arcs = [1 5; 1 13; 1 34; 2 4; 2 5; 2 17; 2 18; 2 20; 2 22; 3 7; 3 9; 3 29; 4 2; 4 32; 5 26; 5 29; 6 16; 7 19; 7 26; 7 31; 7 32; 7 35; 8 22; 8 24; 9 8; 9 11; 9 13; 9 19; 9 30; 10 13; 10 23; 10 28; 11 10; 12 3; 12 13; 12 14; 12 17; 13 22; 13 34; 14 2; 14 27; 15 8; 15 28; 15 30; 15 31; 16 5; 16 6; 16 10; 17 20; 18 8; 18 11; 18 21; 19 3; 19 22; 20 8; 20 11; 20 15; 20 31; 21 4; 21 9; 21 13; 21 28; 21 30; 21 35; 22 3; 22 25; 22 31; 22 32; 23 3; 23 21; 23 24; 23 33; 23 35; 24 6; 25 12; 25 15; 25 32; 25 35; 26 7; 26 8; 26 22; 26 23; 26 30; 26 31; 27 31; 28 4; 28 16; 29 8; 29 15; 29 19; 29 21; 30 6; 30 21; 30 22; 31 10; 31 19; 31 20; 31 32; 32 4; 32 23; 32 26; 32 30; 33 7; 33 31; 34 9; 34 22; 34 23; 34 31]
global d_x = [5.0, 5.0, 4.0, 2.0, 7.0, 1.0, 10.0, 2.0, 7.0, 3.0, 9.0, 10.0, 5.0, 6.0, 9.0, 3.0, 2.0, 4.0, 4.0, 2.0, 3.0, 3.0, 3.0, 9.0, 2.0, 2.0, 1.0, 8.0, 9.0, 7.0, 8.0, 1.0, 6.0, 4.0, 2.0, 6.0, 8.0, 4.0, 4.0, 7.0, 5.0, 10.0, 6.0, 4.0, 2.0, 6.0, 3.0, 1.0, 10.0, 4.0, 2.0, 8.0, 4.0, 4.0, 4.0, 10.0, 3.0, 2.0, 4.0, 5.0, 9.0, 5.0, 9.0, 1.0, 5.0, 6.0, 10.0, 2.0, 2.0, 2.0, 1.0, 8.0, 10.0, 4.0, 10.0, 1.0, 1.0, 7.0, 7.0, 3.0, 6.0, 1.0, 8.0, 9.0, 4.0, 8.0, 8.0, 9.0, 8.0, 4.0, 8.0, 6.0, 10.0, 10.0, 10.0, 9.0, 9.0, 2.0, 9.0, 10.0, 5.0, 1.0, 6.0, 2.0, 4.0, 8.0, 7.0, 7.0]
global b_x = 5
global d_y = [8.0, 9.0, 10.0, 9.0, 10.0, 10.0, 9.0, 2.0, 9.0, 10.0, 4.0, 8.0, 9.0, 9.0, 4.0, 1.0, 8.0, 6.0, 6.0, 8.0, 8.0, 4.0, 9.0, 3.0, 7.0, 4.0, 5.0, 6.0, 2.0, 4.0, 2.0, 3.0, 1.0, 1.0, 7.0, 8.0, 10.0, 3.0, 3.0, 3.0, 5.0, 2.0, 3.0, 5.0, 6.0, 3.0, 6.0, 4.0, 6.0, 9.0, 3.0, 4.0, 7.0, 7.0, 9.0, 9.0, 9.0, 5.0, 9.0, 8.0, 4.0, 2.0, 10.0, 6.0, 6.0, 6.0, 1.0, 8.0, 7.0, 6.0, 5.0, 1.0, 7.0, 5.0, 10.0, 2.0, 7.0, 6.0, 9.0, 3.0, 8.0, 4.0, 9.0, 9.0, 4.0, 10.0, 9.0, 7.0, 6.0, 10.0, 5.0, 9.0, 2.0, 7.0, 1.0, 10.0, 6.0, 8.0, 10.0, 6.0, 7.0, 10.0, 1.0, 6.0, 6.0, 8.0, 2.0, 4.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.994, 0.912, 0.911, 0.953, 0.913, 0.843, 0.99, 0.88, 0.891, 0.921, 0.901, 0.862, 0.808, 0.928, 0.938, 0.898, 0.903, 0.944, 0.918, 0.937, 0.995, 0.927, 0.957, 0.904, 0.846, 0.974, 0.908, 0.854, 0.842, 0.961, 0.848, 0.949, 0.943, 0.83, 0.836, 0.928, 0.946, 0.885, 0.828, 0.95, 0.855, 0.924, 0.991, 0.928, 0.938, 0.825, 0.894, 0.866, 0.957, 0.903, 0.978, 0.977, 0.901, 0.847, 0.846, 0.905, 0.913, 0.823, 0.917, 0.914, 0.847, 0.932, 0.932, 0.924, 0.977, 0.884, 0.973, 0.919, 0.929, 0.949, 0.979, 0.858, 0.923, 0.923, 0.921, 0.841, 0.877, 0.871, 0.915, 0.809, 0.895, 0.934, 0.892, 0.809, 0.95, 0.946, 0.871, 0.809, 0.834, 0.973, 0.92, 0.985, 0.953, 0.941, 0.952, 0.88, 0.801, 0.899, 0.967, 0.851, 0.902, 0.824, 0.96, 0.976, 0.883, 0.85, 0.965, 0.954]
global origin = 1
global destination = 35