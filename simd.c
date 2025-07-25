// CPUでの性能向上 --- マルチコア向上 (OpenMP) + SIMD 化

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include "pcg_variants.h"

#define NX 32                       // 格子のXサイズ
#define NY 32                       // 格子のYサイズ
#define NL 16                       // SIMD化のブロックサイズ
#define N_REPLICAS 16               // レプリカの数
#define N_SAMPLING 1000             // サンプリング数
#define SAMPLE_INTERVAL_SWEEPS 100  // サンプリング間隔
#define THERMALIZATION_SWEEPS 10000 // 初期緩和のスイープ数

// 相互作用の強さ
#define J (1.0f)

/**
 * @brief レプリカの状態を保持する構造体
 */
typedef struct
{
  float *spin;  // スピンの角度θを格納する配列 (サイズ: N_SITES)
  float beta;   // このレプリカの逆温度 β = 1/T
  float energy; // このレプリカの現在のエネルギー
  int nx;       // 格子のX方向のサイズ
  int ny;       // 格子のY方向のサイズ
} Replica;

/**
 * @brief 系全体のエネルギーを計算する
 * @param rep 計算対象のレプリカへのポインタ
 * @return 計算された全エネルギー
 * @note 周期境界条件を適用します
 */
float calculate_total_energy(Replica *rep)
{
  float total_energy = 0.0f;
  const int N_SITES = rep->nx * rep->ny;

  for (int i = 0; i < N_SITES; ++i)
  {
    int ix = i % rep->nx;
    int iy = i / rep->nx;

    // 右の隣接スピンとの相互作用
    int right_neighbor_idx = iy * rep->nx + (ix + 1) % rep->nx;
    total_energy -= J * cosf(rep->spin[i] - rep->spin[right_neighbor_idx]);

    // 下の隣接スピンとの相互作用
    int down_neighbor_idx = ((iy + 1) % rep->ny) * rep->nx + ix;
    total_energy -= J * cosf(rep->spin[i] - rep->spin[down_neighbor_idx]);
  }
  return total_energy;
}

/**
 * @brief パラレルテンパレリング用のレプリカ群を初期化する
 * @param replicas レプリカの配列へのポインタ
 * @param n_replicas レプリカの数
 * @param nx 格子のX方向のサイズ
 * @param ny 格子のY方向のサイズ
 * @param betas 各レプリカに設定する逆温度の配列
 * @note 乱数シードは事前に設定しておく必要があります
 */
void initialize_replicas(Replica *replicas, int n_replicas, int nx, int ny, const float *betas)
{
  const int N_SITES = nx * ny;
  for (int i = 0; i < n_replicas; ++i)
  {
    replicas[i].spin = (float *)malloc(sizeof(float) * N_SITES);
    if (replicas[i].spin == NULL)
    {
      fprintf(stderr, "Error: Failed to allocate memory for spins.\n");
      exit(EXIT_FAILURE);
    }

    replicas[i].beta = betas[i];
    replicas[i].nx = nx;
    replicas[i].ny = ny;

    // スピンをランダムな角度 [0, 2π) で初期化
    for (int j = 0; j < N_SITES; ++j)
    {
      replicas[i].spin[j] = ((float)rand() / (float)RAND_MAX) * 2.0f * M_PI;
    }

    // 初期エネルギーを計算
    replicas[i].energy = calculate_total_energy(&replicas[i]);
  }
}

/**
 * @brief 1つのレプリカに対して1モンテカルロスイープを実行する (メトロポリス法)
 * @param rep 更新するレプリカへのポインタ
 * @note 全てのスピンサイトを更新します
 */
void metropolis_sweep(Replica *rep, pcg32_random_t *rng)
{
  const int N_SITES = rep->nx * rep->ny;
  const float DELTA = 1.0f; // スピン角度の更新幅

  int odd, site_idx, idx, ix, iy, r; // loop variables

  // block variables for SIMD
  float spin[NL], new_spin[NL];                                         // spin block
  float neighborsR[NL], neighborsL[NL], neighborsD[NL], neighborsU[NL]; // neighbors blocks
  float delta_e[NL];                                                    // energy change and acceptance probability
  float rand_val_acc[NL], rand_val_delta[NL];                           // random values for acceptance and delta
  int accept[NL];                                                       // acceptance flags

  // update each site in blocks of NL
  for (odd = 0; odd < 2; odd++)
  {
    // odd = 0 : even sites
    // odd = 1 : odd sites
    for (site_idx = odd; site_idx < N_SITES; site_idx += NL * 2)
    {
      // gather load
      for (r = 0; r < NL; ++r)
      {
        idx = site_idx + r * 2;
        ix = idx % rep->nx;
        iy = idx / rep->nx;

        // load spin
        spin[r] = rep->spin[idx];

        // load neighbors
        neighborsR[r] = rep->spin[iy * rep->nx + (ix + 1) % rep->nx];             // Right
        neighborsL[r] = rep->spin[iy * rep->nx + (ix - 1 + rep->nx) % rep->nx];   // Left
        neighborsD[r] = rep->spin[((iy + 1) % rep->ny) * rep->nx + ix];           // Down
        neighborsU[r] = rep->spin[((iy - 1 + rep->ny) % rep->ny) * rep->nx + ix]; // Up

        rand_val_delta[r] = (float)pcg32_random_r(rng) / (float)UINT32_MAX;
        rand_val_acc[r] = (float)pcg32_random_r(rng) / (float)UINT32_MAX;
      }

#pragma omp simd
      for (int r = 0; r < NL; ++r)
      {
        new_spin[r] = fmodf(spin[r] + (rand_val_delta[r] - 0.5f) * DELTA + 2.0f * M_PI, 2.0f * M_PI);

        // calculate energy change
        delta_e[r] = 0.0f;
        // calculate energy change with neighbors
        delta_e[r] += -J * (cosf(new_spin[r] - neighborsR[r]) - cosf(spin[r] - neighborsR[r]));
        delta_e[r] += -J * (cosf(new_spin[r] - neighborsL[r]) - cosf(spin[r] - neighborsL[r]));
        delta_e[r] += -J * (cosf(new_spin[r] - neighborsD[r]) - cosf(spin[r] - neighborsD[r]));
        delta_e[r] += -J * (cosf(new_spin[r] - neighborsU[r]) - cosf(spin[r] - neighborsU[r]));

        accept[r] = (delta_e[r] < 0.0f || rand_val_acc[r] < expf(-rep->beta * delta_e[r]));
      }

      // scatter store
      for (int r = 0; r < NL; ++r)
      {
        int idx = site_idx + r * 2;
        if (accept[r])
        {
          // accept the new spin
          rep->spin[idx] = new_spin[r];
          rep->energy += delta_e[r];
        }
      }
    }
  }
}

/**
 * @brief レプリカ交換を実行する
 * @param replicas レプリカの配列
 * @param n_replicas レプリカの数
 */
void replica_exchange(Replica *replicas, int n_replicas)
{
  // 隣接するレプリカのペア (i, i+1) をランダムに選択
  int i = (rand() % (n_replicas - 1));

  float delta_beta = replicas[i].beta - replicas[i + 1].beta;
  float delta_energy = replicas[i].energy - replicas[i + 1].energy;

  // 交換確率を計算: min(1, exp(Δβ * ΔE))
  float acceptance_prob = expf(delta_beta * delta_energy);

  if (rand() < (float)RAND_MAX * acceptance_prob)
  {
    // スピン配列とエネルギーを交換
    float *temp_spin = replicas[i].spin;
    replicas[i].spin = replicas[i + 1].spin;
    replicas[i + 1].spin = temp_spin;

    float temp_energy = replicas[i].energy;
    replicas[i].energy = replicas[i + 1].energy;
    replicas[i + 1].energy = temp_energy;
  }
}

/**
 * @brief 確保したメモリを解放する
 * @param replicas レプリカの配列
 * @param n_replicas レプリカの数
 */
void free_replicas(Replica *replicas, int n_replicas)
{
  for (int i = 0; i < n_replicas; ++i)
  {
    if (replicas[i].spin != NULL)
    {
      free(replicas[i].spin);
      replicas[i].spin = NULL;
    }
  }
}

int main(int argc, char *argv[])
{
  // 1. パラメータ設定
  const int N_SITES = NX * NY;
  Replica replicas[N_REPLICAS];
  float betas[N_REPLICAS];
  int i, sweep;

  float min_beta = 0.4f;
  float max_beta = 1.9f;

  struct timeval start_time, end_time;
  double wall_time;

  FILE *fp;

  // start time measurement
  gettimeofday(&start_time, NULL);

  if (N_SITES % (NL * 2) != 0)
  {
    fprintf(stderr, "Error: N_SITES must be a multiple of %d.\n", NL * 2);
    return EXIT_FAILURE;
  }

  if (argc != 2)
  {
    fprintf(stderr, "Usage: %s <output_file>\n", argv[0]);
    return EXIT_FAILURE;
  }

  // 出力ファイルを開く
  fp = fopen(argv[1], "w");
  if (fp == NULL)
  {
    fprintf(stderr, "Error: Could not open output file %s\n", argv[1]);
    return EXIT_FAILURE;
  }

  fprintf(fp, "<version\nsimd\nversion>\n");

  fprintf(fp, "<input_parameters\n");
  fprintf(fp, "N_SITES=%d N_REPLICAS=%d NX=%d NY=%d N_SAMPLING=%d SAMPLE_INTERVAL_SWEEPS=%d THERMALIZATION_SWEEPS=%d\n",
          N_SITES, N_REPLICAS, NX, NY, N_SAMPLING, SAMPLE_INTERVAL_SWEEPS, THERMALIZATION_SWEEPS);
  fprintf(fp, "min_beta=%.6f max_beta=%.6f\n", min_beta, max_beta);
  fprintf(fp, "input_parameters>\n");

  for (i = 0; i < N_REPLICAS; ++i)
  {
    betas[i] = min_beta + (max_beta - min_beta) * i / (N_REPLICAS - 1);
  }

  // 2. 初期化
  srand(48);
  int max_threads = omp_get_max_threads();
  pcg32_random_t rngs[max_threads];
  for (i = 0; i < max_threads; ++i)
  {
    pcg32_srandom_r(&rngs[i], time(NULL), (intptr_t)&rngs[i]);
  }
  initialize_replicas(replicas, N_REPLICAS, NX, NY, betas);

  // 全レプリカの物理量を格納する配列
  double total_Ex[N_REPLICAS] = {0.0};
  double total_Jx_squared[N_REPLICAS] = {0.0};
  double total_Ey[N_REPLICAS] = {0.0};
  double total_Jy_squared[N_REPLICAS] = {0.0};
  long measurement_count = 0;

  // 3. シミュレーションループ
  printf("<simulation_progress\n");
  sweep = 0;
  while (measurement_count < N_SAMPLING)
  {
    // 各レプリカでメトロポリス更新
#pragma omp parallel for
    for (i = 0; i < N_REPLICAS; ++i)
    {
      metropolis_sweep(&replicas[i], &rngs[omp_get_thread_num()]);
    }

    // レプリカ交換
    for (i = 0; i < N_REPLICAS; ++i)
    { // 交換頻度を上げるためにループ
      replica_exchange(replicas, N_REPLICAS);
    }

    // 物理量の測定(初期緩和後)
    if (sweep > THERMALIZATION_SWEEPS && (sweep - THERMALIZATION_SWEEPS) % SAMPLE_INTERVAL_SWEEPS == 0)
    {
#pragma omp parallel for
      for (i = 0; i < N_REPLICAS; ++i)
      {
        float E_x_sweep = 0.0f;
        float J_x_sweep = 0.0f;
        float E_y_sweep = 0.0f;
        float J_y_sweep = 0.0f;

#pragma omp simd reduction(+ : E_x_sweep, J_x_sweep, E_y_sweep, J_y_sweep)
        for (int j = 0; j < N_SITES; ++j)
        {
          int ix = j % NX;
          int iy = j / NX;

          // x方向 (右の隣人)
          int right_neighbor_idx = iy * NX + (ix + 1) % NX;
          float delta_theta_x = replicas[i].spin[j] - replicas[i].spin[right_neighbor_idx];
          E_x_sweep += J * cosf(delta_theta_x);
          J_x_sweep += J * sinf(delta_theta_x);

          // y方向 (下の隣人)
          int down_neighbor_idx = ((iy + 1) % NY) * NX + ix;
          float delta_theta_y = replicas[i].spin[j] - replicas[i].spin[down_neighbor_idx];
          E_y_sweep += J * cosf(delta_theta_y);
          J_y_sweep += J * sinf(delta_theta_y);
        }
        total_Ex[i] += E_x_sweep;
        total_Jx_squared[i] += (double)J_x_sweep * J_x_sweep;
        total_Ey[i] += E_y_sweep;
        total_Jy_squared[i] += (double)J_y_sweep * J_y_sweep;
      }
      measurement_count++;
      if (measurement_count % (N_SAMPLING / 10) == 0)
      {
        printf("sweep=%d measurement_count=%ld total_Ex=%f total_Jx_squared=%f total_Ey=%f total_Jy_squared=%f\n",
               sweep, measurement_count, total_Ex[0], total_Jx_squared[0], total_Ey[0], total_Jy_squared[0]);
      }
    }

    sweep++;
  }
  printf("simulation_progress>\n");

  // 最終結果の計算と表示
  fprintf(fp, "<helicity_modulus_results\n");
  for (i = 0; i < N_REPLICAS; ++i)
  {
    if (measurement_count > 0)
    {
      double avg_Ex = total_Ex[i] / measurement_count;
      double avg_Jx_squared = total_Jx_squared[i] / measurement_count;
      double upsilon_x = (avg_Ex / N_SITES) - (replicas[i].beta / N_SITES) * avg_Jx_squared;

      double avg_Ey = total_Ey[i] / measurement_count;
      double avg_Jy_squared = total_Jy_squared[i] / measurement_count;
      double upsilon_y = (avg_Ey / N_SITES) - (replicas[i].beta / N_SITES) * avg_Jy_squared;

      double helicity_modulus = (upsilon_x + upsilon_y) / 2.0;
      double energy_per_site = replicas[i].energy / (double)N_SITES;
      fprintf(fp, "replica=%d beta=%.6f energy_per_site=%.6f helicity_modulus=%.6f\n",
              i, replicas[i].beta, energy_per_site, helicity_modulus);
    }
  }
  fprintf(fp, "helicity_modulus_results>\n");

  // 4. メモリ解放
  free_replicas(replicas, N_REPLICAS);

  // end time measurement
  gettimeofday(&end_time, NULL);
  wall_time = (end_time.tv_sec - start_time.tv_sec) +
              (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
  fprintf(fp, "<execution_time\n%.4f\nexecution_time>\n", wall_time);

  fclose(fp);

  return 0;
}
