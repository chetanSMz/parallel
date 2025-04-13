#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>

double getRandomNumber(double min, double max) {
    double scale = rand() / (double)RAND_MAX;
    double val = min + scale * (max - min);
    return round(val * 10000.0) / 10000.0;
}

int main(int argc, char *argv[]) {
    srand((unsigned int)time(NULL));

    int row, col;
    row=atoi(argv[1]);
    col=atoi(argv[2]);

    int nx, ny, nz;
    nx=atoi(argv[3]);
    ny=atoi(argv[4]);
    nz=atoi(argv[5]);

    float **v = (float **)malloc(col * sizeof(float *));
    double **fa = (double **)malloc(col * sizeof(double *));
    for (int i = 0; i < col; i++) {
        v[i] = (float *)malloc(row * sizeof(float));
        fa[i] = (double *)malloc(4 * sizeof(double));
    }

    int prev_col = col;
    int curr_idx = 0;

    int dx[] = {1, 0, -1, 0, 0, 0};
    int dy[] = {0, 1, 0, -1, 0, 0};
    int dz[] = {0, 0, 0, 0, 1, -1};

    while (col--) {
        float *tmp = (float *)malloc(row * sizeof(float));
        double ***dp = (double ***)malloc(nx * sizeof(double **));
        for (int x = 0; x < nx; x++) {
            dp[x] = (double **)malloc(ny * sizeof(double *));
            for (int y = 0; y < ny; y++) {
                dp[x][y] = (double *)malloc(nz * sizeof(double));
            }
        }

        double mx1 = -DBL_MAX, mi1 = DBL_MAX;
        int count = 0;

        for (int z = 0; z < nz; z++) {
            for (int y = 0; y < ny; y++) {
                for (int x = 0; x < nx; x++) {
                    double num = getRandomNumber(-100.0, 100.0);
                    dp[x][y][z] = num;
                    if (num > mx1) mx1 = num;
                    if (num < mi1) mi1 = num;
                    tmp[count++] = (float)num;
                }
            }
        }

        int lmax = 0, lmi = 0;

        for (int z = 0; z < nz; z++) {
            for (int y = 0; y < ny; y++) {
                for (int x = 0; x < nx; x++) {
                    double curr = dp[x][y][z];
                    double mx = curr, mi = curr;int c1=1;int c2=1;

                    for (int i = 0; i < 6; i++) {
                        int x_new = x + dx[i];
                        int y_new = y + dy[i];
                        int z_new = z + dz[i];

                        if (x_new >= 0 && x_new < nx && y_new >= 0 && y_new < ny && z_new >= 0 && z_new < nz) {
                            double val = dp[x_new][y_new][z_new];
                            if(val>mx){
                                mx=val;c1=1;
                            }
                            else if(val==mx)c1++;
                            if(val<mi){
                                mi=val;c2=1;
                            }
                            else if(val==mi){
                                c2++;
                            }
                        }
                    }

                    if (curr == mx && c1==1) lmax++;
                    if (curr == mi && c2==1) lmi++;
                }
            }
        }

        fa[curr_idx][0] = (double)lmi;
        fa[curr_idx][1] = (double)lmax;
        fa[curr_idx][2] = mx1;
        fa[curr_idx][3] = mi1;

        for (int i = 0; i < row; i++) {
            v[curr_idx][i] = tmp[i];
        }

        // Cleanup
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                free(dp[x][y]);
            }
            free(dp[x]);
        }
        free(dp);
        free(tmp);

        curr_idx++;
    }

    FILE *mybin = fopen("input.bin", "wb");
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < prev_col; j++) {
            float val = v[j][i];
            fwrite(&val, sizeof(float), 1, mybin);
        }
    }
    fclose(mybin);

    FILE *myfile2 = fopen("output.txt", "w");
    for (int i = 0; i < prev_col; i++) {
        fprintf(myfile2, "%d\n", i);
        for (int j = 0; j < 4; j++) {
            fprintf(myfile2, "%.4f ", fa[i][j]);
        }
        fprintf(myfile2, "\n\n");
    }
    fclose(myfile2);

    // Cleanup
    for (int i = 0; i < prev_col; i++) {
        free(v[i]);
        free(fa[i]);
    }
    free(v);
    free(fa);

    return 0;
}
