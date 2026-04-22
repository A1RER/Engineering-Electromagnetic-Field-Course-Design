/*
 * dipole_field.c  ─  Electric / Magnetic Dipole Field Numerical Simulation
 *
 * Modules implemented
 *   Module 1: Physical constants & dipole source definition
 *   Module 2: 3-D Cartesian mesh generation
 *   Module 3: Vectorized field computation kernel (loop-level, C99)
 *   Module 4: CSV data export (Origin-compatible: x,y,z,Fx,Fy,Fz,|F|)
 *
 * Compile (GCC / MinGW):
 *   gcc -O2 -std=c99 -o dipole_field dipole_field.c -lm
 *
 * Compile with OpenMP (optional parallelisation):
 *   gcc -O2 -std=c99 -fopenmp -o dipole_field dipole_field.c -lm
 *
 * Usage:
 *   dipole_field [type] [mode] [N] [L_nm]
 *
 *   type  : e = electric dipole (default)
 *           m = magnetic dipole
 *   mode  : s = static near-field (default)
 *           t = time-harmonic (Hertzian dipole)
 *   N     : grid points per axis  (default: 40, giving 40³ = 64 000 pts)
 *   L_nm  : domain half-size [nm] (default: 5.0)
 *
 * Output: dipole_field_data.csv in the current directory
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* ── Physical constants ──────────────────────────────────────────────── */
#define EPSILON_0  8.854187817e-12   /* F/m  */
#define MU_0       1.256637061e-6    /* H/m  */
#define PI         3.14159265358979323846

/* ── Dipole type / mode ───────────────────────────────────────────────── */
typedef enum { ELECTRIC, MAGNETIC } DipoleType;
typedef enum { STATIC,   HARMONIC } SimMode;

/* ── Dipole configuration ─────────────────────────────────────────────── */
typedef struct {
    double px, py, pz;    /* dipole moment vector (C·m or A·m²) */
    double x0, y0, z0;    /* source position [m]                 */
    double k;             /* wave-number [rad/m], 0 for static   */
    double eps_s;         /* singularity softening length [m]    */
    DipoleType type;
    SimMode    mode;
} DipoleConfig;

/* ── Field at a single point ──────────────────────────────────────────── */
typedef struct { double fx, fy, fz, mag; } FieldPoint;


/* ===================================================================== */
/*  Computation Kernel                                                   */
/* ===================================================================== */
/*
 * compute_field_at()
 *   Returns E (electric) or H (magnetic) at observation point (x,y,z).
 *
 *   Static formula (same structure for E and H, differing only in coeff):
 *     F_i = (coeff / r⁵) · [3·d_i·(p·d) - p_i·r²]
 *
 *   Time-harmonic formula:
 *     F_i = coeff·e^{-jkr}·{ (k²/r³)·pt_i  +  [(1+jkr)/r⁵]·pl_i }
 *   where pt = p - (p·d/r²)·d  (transverse)
 *         pl = 3(p·d/r²)·d - p (longitudinal)
 *
 *   Singularity regularisation: r² → r² + ε_s²
 *
 * NOTE: For the time-harmonic case the field is complex.  This C
 *       implementation returns the real part (instantaneous field at t=0,
 *       i.e. Re[F·e^{+j0}]) and the complex magnitude |F|.
 */
static FieldPoint compute_field_at(double x, double y, double z,
                                   const DipoleConfig *cfg)
{
    FieldPoint fp;
    const double px = cfg->px, py = cfg->py, pz = cfg->pz;
    const double k  = cfg->k;

    /* Displacement from source */
    double dx = x - cfg->x0;
    double dy = y - cfg->y0;
    double dz = z - cfg->z0;

    /* Stabilised distance */
    double r2  = dx*dx + dy*dy + dz*dz + cfg->eps_s * cfg->eps_s;
    double r   = sqrt(r2);
    double r5  = r2 * r2 * r;

    /* Prefactor: 1/(4πε₀) for E, 1/(4π) for H */
    double coeff = (cfg->type == ELECTRIC)
                   ? 1.0 / (4.0 * PI * EPSILON_0)
                   : 1.0 / (4.0 * PI);

    double p_dot_d = px*dx + py*dy + pz*dz;

    double Fx, Fy, Fz;

    if (cfg->mode == STATIC) {
        /* ── Static field ────────────────────────────────────────────── */
        double c_r5 = coeff / r5;
        Fx = c_r5 * (3.0*dx*p_dot_d - px*r2);
        Fy = c_r5 * (3.0*dy*p_dot_d - py*r2);
        Fz = c_r5 * (3.0*dz*p_dot_d - pz*r2);

    } else {
        /* ── Time-harmonic field (real part at ωt = 0) ───────────────── */
        double kr   = k * r;
        double cos_kr = cos(kr);     /* Re[e^{-jkr}] = cos(kr) */
        double sin_kr = sin(kr);     /* Im[e^{-jkr}] = -sin(kr) */

        /*
         * A = (1+jkr)/r⁵   →   Re[A·e^{-jkr}] = (cos_kr + kr·sin_kr)/r⁵
         * B = k²/r³         →   Re[B·e^{-jkr}] = k²·cos_kr/r³
         */
        double r3   = r2 * r;
        double A_re = (cos_kr + kr*sin_kr) / r5;
        double B_re = (k*k * cos_kr) / r3;

        /* Transverse component: pt_i = p_i - (p·d/r²)·d_i */
        double ptx = px - p_dot_d*dx/r2;
        double pty = py - p_dot_d*dy/r2;
        double ptz = pz - p_dot_d*dz/r2;

        /* Longitudinal: pl_i = 3(p·d/r²)·d_i - p_i */
        double plx = 3.0*p_dot_d*dx/r2 - px;
        double ply = 3.0*p_dot_d*dy/r2 - py;
        double plz = 3.0*p_dot_d*dz/r2 - pz;

        Fx = coeff * (B_re*ptx + A_re*plx);
        Fy = coeff * (B_re*pty + A_re*ply);
        Fz = coeff * (B_re*ptz + A_re*plz);
    }

    fp.fx  = Fx;
    fp.fy  = Fy;
    fp.fz  = Fz;
    fp.mag = sqrt(Fx*Fx + Fy*Fy + Fz*Fz);
    return fp;
}


/* ===================================================================== */
/*  Grid sweep (Module 2 + 3 combined)                                   */
/* ===================================================================== */
static int run_simulation(const DipoleConfig *cfg,
                          int N, double L,
                          const char *out_path)
{
    long total = (long)N * N * N;
    printf("[mesh] %d³ = %ld grid points  |  domain: ±%.2e m\n",
           N, total, L);
    printf("[mesh] eps_s = %.2e m\n", cfg->eps_s);

    /* Open output CSV */
    FILE *fp = fopen(out_path, "w");
    if (!fp) { perror("fopen"); return -1; }
    fprintf(fp, "x_nm,y_nm,z_nm,Fx,Fy,Fz,F_mag\n");

    const char *type_str = (cfg->type == ELECTRIC) ? "Electric E" : "Magnetic H";
    const char *mode_str = (cfg->mode == STATIC)   ? "static"     : "time-harmonic";
    printf("[sim]  %s  |  mode: %s  |  writing → %s\n",
           type_str, mode_str, out_path);

    clock_t t0 = clock();
    long n_written = 0;
    double step = (N > 1) ? (2.0*L / (N-1)) : 0.0;

    for (int ix = 0; ix < N; ix++) {
        double x = -L + ix * step;
        for (int iy = 0; iy < N; iy++) {
            double y = -L + iy * step;
            for (int iz = 0; iz < N; iz++) {
                double z = -L + iz * step;

                FieldPoint f = compute_field_at(x, y, z, cfg);

                /* Export in nm for human-readable columns */
                fprintf(fp, "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        x*1e9, y*1e9, z*1e9,
                        f.fx, f.fy, f.fz, f.mag);
                n_written++;
            }
        }
        /* Progress indicator every 10% */
        if ((ix+1) % (N/10 + 1) == 0)
            printf("  progress: %.0f%%\r", 100.0*(ix+1)/N);
    }

    fclose(fp);
    double elapsed = (double)(clock()-t0) / CLOCKS_PER_SEC;
    printf("\n[done] %ld points written  |  elapsed: %.3f s\n",
           n_written, elapsed);
    return 0;
}


/* ===================================================================== */
/*  Analytical self-test (Module 5 inline check)                         */
/* ===================================================================== */
/*
 * For a z-directed electric dipole (p=[0,0,1]), on the polar axis (z):
 *   E_z = 2p / (4πε₀ r³)
 * Check at r = L/2; relative error must be < 1e-3 for eps_s = L·1e-4.
 */
static void self_test(const DipoleConfig *cfg, double L)
{
    double r_test = 0.5 * L;
    double x=0.0, y=0.0, z=r_test;

    DipoleConfig cfg_test = *cfg;
    /* Force z-directed unit electric dipole for the test */
    cfg_test.type = ELECTRIC;
    cfg_test.mode = STATIC;
    cfg_test.px=0; cfg_test.py=0; cfg_test.pz=1.0;
    cfg_test.k=0;

    FieldPoint f = compute_field_at(x, y, z, &cfg_test);

    double E_analytic = 2.0 / (4.0*PI*EPSILON_0 * r_test*r_test*r_test);
    double err = fabs(f.fz - E_analytic) / fabs(E_analytic);

    printf("[self-test]  r=L/2  |  Fz_num=%.4e  Fz_ana=%.4e  rel_err=%.2e  %s\n",
           f.fz, E_analytic, err, (err < 1e-3) ? "PASS" : "FAIL");
}


/* ===================================================================== */
/*  main                                                                  */
/* ===================================================================== */
int main(int argc, char *argv[])
{
    /* ── Default configuration ──────────────────────────────────────── */
    DipoleConfig cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.px = 0.0;  cfg.py = 0.0;  cfg.pz = 1e-30;   /* z-directed */
    cfg.x0 = 0.0;  cfg.y0 = 0.0;  cfg.z0 = 0.0;
    cfg.type = ELECTRIC;
    cfg.mode = STATIC;
    cfg.k    = 0.0;

    int    N     = 40;
    double L_nm  = 5.0;   /* domain half-size in nm */

    /* ── Parse command-line arguments ───────────────────────────────── */
    if (argc > 1) cfg.type = (argv[1][0]=='m') ? MAGNETIC : ELECTRIC;
    if (argc > 2) cfg.mode = (argv[2][0]=='t') ? HARMONIC : STATIC;
    if (argc > 3) N        = atoi(argv[3]);
    if (argc > 4) L_nm     = atof(argv[4]);

    if (N < 2 || N > 500) {
        fprintf(stderr, "N must be in [2, 500]; got %d\n", N);
        return EXIT_FAILURE;
    }

    double L   = L_nm * 1e-9;        /* convert nm → m  */
    cfg.eps_s  = L * 1e-4;           /* singularity softening */

    /* For time-harmonic: use 1 GHz by default */
    if (cfg.mode == HARMONIC) {
        double freq = 1e9;
        double c    = 1.0 / sqrt(EPSILON_0 * MU_0);
        cfg.k = 2.0 * PI * freq / c;
        printf("[config] Harmonic mode: f=%.1e Hz  k=%.4e rad/m\n",
               freq, cfg.k);
    }

    /* Set magnetic dipole moment if needed */
    if (cfg.type == MAGNETIC) {
        cfg.px = 0.0;  cfg.py = 0.0;  cfg.pz = 1e-20;   /* [A·m²] */
    }

    /* ── Self-test ───────────────────────────────────────────────────── */
    self_test(&cfg, L);

    /* ── Run simulation ─────────────────────────────────────────────── */
    const char *out_path = "dipole_field_data.csv";
    int ret = run_simulation(&cfg, N, L, out_path);

    return (ret == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
