global arcs = [1 5; 1 11; 1 23; 1 25; 1 29; 1 35; 2 7; 2 9; 3 9; 3 23; 3 27; 3 29; 3 39; 4 26; 4 35; 5 11; 5 18; 5 32; 5 38; 6 14; 6 17; 7 22; 7 37; 8 2; 8 6; 8 9; 8 27; 8 29; 8 33; 9 27; 10 9; 10 12; 10 17; 10 21; 10 32; 10 35; 10 37; 11 7; 11 13; 11 21; 11 26; 12 22; 12 28; 12 30; 12 36; 12 38; 13 8; 13 31; 14 5; 14 8; 14 9; 14 12; 14 18; 14 21; 14 40; 15 16; 15 17; 15 24; 15 26; 15 39; 16 3; 16 8; 16 14; 16 29; 16 33; 16 38; 17 8; 17 33; 17 35; 17 37; 18 3; 18 13; 18 17; 18 36; 19 7; 19 20; 19 21; 19 28; 20 4; 20 22; 20 24; 21 10; 21 23; 21 24; 21 29; 21 32; 22 5; 22 8; 22 14; 22 19; 22 20; 22 28; 22 31; 22 34; 22 37; 22 40; 23 11; 23 28; 23 30; 24 5; 24 14; 24 23; 24 36; 25 8; 25 31; 25 35; 25 36; 26 17; 26 23; 27 6; 27 16; 27 22; 27 39; 28 3; 28 4; 28 7; 28 12; 28 21; 28 23; 28 26; 29 4; 29 13; 29 19; 30 22; 31 21; 31 23; 31 32; 31 37; 31 40; 32 4; 32 20; 32 22; 32 25; 32 31; 33 26; 34 2; 34 5; 34 15; 34 31; 34 36; 35 15; 35 22; 35 33; 36 3; 36 24; 37 10; 37 11; 37 12; 37 17; 37 20; 37 29; 37 32; 38 5; 38 17; 38 19; 38 28; 38 32; 38 33; 38 37; 39 23; 39 27; 39 32]
global d_x = [1.0, 7.0, 7.0, 6.0, 7.0, 2.0, 8.0, 9.0, 6.0, 4.0, 3.0, 10.0, 4.0, 7.0, 9.0, 8.0, 3.0, 6.0, 9.0, 4.0, 7.0, 3.0, 6.0, 10.0, 3.0, 1.0, 6.0, 8.0, 2.0, 6.0, 6.0, 4.0, 4.0, 6.0, 1.0, 7.0, 9.0, 4.0, 10.0, 7.0, 6.0, 10.0, 7.0, 4.0, 3.0, 8.0, 10.0, 9.0, 8.0, 5.0, 7.0, 2.0, 3.0, 9.0, 8.0, 6.0, 7.0, 7.0, 5.0, 10.0, 9.0, 8.0, 1.0, 4.0, 10.0, 3.0, 7.0, 2.0, 3.0, 8.0, 2.0, 4.0, 3.0, 3.0, 10.0, 1.0, 2.0, 3.0, 1.0, 10.0, 2.0, 10.0, 1.0, 5.0, 8.0, 10.0, 9.0, 4.0, 8.0, 3.0, 2.0, 10.0, 5.0, 10.0, 7.0, 3.0, 5.0, 8.0, 2.0, 5.0, 7.0, 8.0, 1.0, 7.0, 4.0, 1.0, 9.0, 1.0, 3.0, 5.0, 5.0, 10.0, 10.0, 3.0, 10.0, 6.0, 4.0, 6.0, 9.0, 7.0, 3.0, 3.0, 10.0, 10.0, 2.0, 4.0, 2.0, 2.0, 2.0, 10.0, 1.0, 9.0, 10.0, 3.0, 10.0, 5.0, 2.0, 9.0, 8.0, 1.0, 9.0, 4.0, 8.0, 3.0, 7.0, 6.0, 6.0, 3.0, 7.0, 5.0, 7.0, 8.0, 7.0, 9.0, 2.0, 1.0, 3.0, 5.0, 5.0, 1.0, 2.0, 1.0]
global b_x = 5
global d_y = [9.0, 2.0, 3.0, 3.0, 6.0, 10.0, 7.0, 1.0, 2.0, 8.0, 10.0, 1.0, 1.0, 5.0, 6.0, 1.0, 9.0, 4.0, 10.0, 7.0, 9.0, 7.0, 1.0, 5.0, 7.0, 2.0, 8.0, 8.0, 5.0, 4.0, 2.0, 9.0, 10.0, 8.0, 8.0, 10.0, 1.0, 7.0, 5.0, 2.0, 8.0, 3.0, 5.0, 6.0, 10.0, 8.0, 4.0, 6.0, 7.0, 10.0, 2.0, 6.0, 5.0, 9.0, 2.0, 8.0, 5.0, 10.0, 10.0, 1.0, 7.0, 5.0, 6.0, 4.0, 8.0, 5.0, 9.0, 9.0, 2.0, 10.0, 1.0, 10.0, 3.0, 7.0, 6.0, 4.0, 2.0, 6.0, 2.0, 3.0, 8.0, 3.0, 5.0, 9.0, 10.0, 3.0, 5.0, 4.0, 10.0, 3.0, 10.0, 7.0, 9.0, 4.0, 10.0, 2.0, 2.0, 3.0, 7.0, 10.0, 7.0, 2.0, 9.0, 3.0, 1.0, 5.0, 3.0, 9.0, 6.0, 5.0, 1.0, 10.0, 9.0, 5.0, 3.0, 1.0, 10.0, 9.0, 3.0, 2.0, 5.0, 8.0, 8.0, 8.0, 3.0, 5.0, 10.0, 4.0, 4.0, 2.0, 3.0, 2.0, 3.0, 7.0, 9.0, 8.0, 7.0, 1.0, 9.0, 1.0, 10.0, 1.0, 3.0, 10.0, 3.0, 6.0, 4.0, 7.0, 10.0, 10.0, 9.0, 7.0, 4.0, 4.0, 2.0, 2.0, 5.0, 9.0, 4.0, 2.0, 1.0, 3.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.864, 0.896, 0.812, 0.879, 0.966, 0.969, 0.956, 0.887, 0.939, 0.852, 0.888, 0.964, 0.809, 0.894, 0.974, 0.834, 0.959, 0.801, 0.844, 0.81, 0.869, 0.806, 0.953, 0.915, 0.826, 0.971, 0.981, 0.928, 0.977, 0.986, 0.945, 0.826, 0.849, 0.881, 0.964, 0.859, 0.82, 0.843, 0.866, 0.927, 0.974, 0.96, 0.88, 0.865, 0.965, 0.994, 0.973, 0.905, 0.994, 0.868, 0.862, 0.96, 0.854, 0.813, 0.826, 0.885, 0.96, 0.873, 0.836, 0.977, 0.865, 0.835, 0.931, 0.868, 0.937, 0.866, 0.929, 0.805, 0.985, 0.967, 0.909, 0.896, 0.909, 0.853, 0.935, 0.979, 0.84, 0.891, 0.808, 0.991, 0.88, 0.963, 0.816, 0.904, 0.94, 0.998, 0.955, 0.881, 0.824, 0.99, 0.853, 0.812, 0.949, 0.803, 0.968, 0.993, 0.872, 0.833, 0.845, 0.998, 0.882, 0.925, 0.923, 0.896, 0.91, 0.841, 0.888, 0.887, 0.916, 0.911, 0.9, 0.886, 0.859, 0.988, 0.895, 0.887, 0.846, 0.863, 0.985, 0.826, 0.909, 0.964, 0.88, 0.8, 0.908, 0.863, 0.916, 0.822, 0.895, 0.808, 0.96, 0.889, 0.933, 0.917, 0.807, 0.946, 0.922, 0.802, 0.931, 0.943, 0.86, 0.878, 0.816, 0.956, 0.819, 0.847, 0.995, 0.858, 0.868, 0.993, 0.898, 0.896, 0.843, 0.9, 0.876, 0.811, 0.961, 0.806, 0.884, 0.888, 0.856, 0.94]
global origin = 1
global destination = 40