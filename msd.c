#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 216
// Ahora bucle en rhos
double rhos[] = {0.2};
int nrho = 5;
#define KB 1.0
#define mass 1.0
#define cutoff 2.5
#define N_EQ 2000
#define nu 0.5 // frec de colision Andersen
#define pi 3.141592653589793
#define NSIM 10000

// Prototipos
static inline double SafeRandom();
void BoxMuller(double* num1, double* num2);

void init_sc(double *x, double *y, double *z, int n_cells, double L);
void init_velocities(double *vx, double *vy, double *vz, double T);
void save_xyz(const char *filename, double *x, double *y, double *z);
void save_velocities(const char *filename, double *vx, double *vy, double *vz);
void compute_LJ_forces_pbc(double *x, double *y, double *z,
                           double *fx, double *fy, double *fz,
                            double L, double *virial);
double compute_LJ_potential(double *x, double *y, double *z,
                            double L);
double kinetic_energy(double *vx, double *vy, double *vz);
double minimum_image(double dx, double L);
void apply_pbc_positions(double *x, double *y, double *z, double L);
void verlet_velocity_step(double *x, double *y, double *z,
                          double *vx, double *vy, double *vz,
                          double *fx, double *fy, double *fz,
                           double L, double dt,double *virial);
void verlet_velocity_step_thermostat(double *x, double *y, double *z,
                          double *vx, double *vy, double *vz,
                          double *fx, double *fy, double *fz,
                          double L, double dt, double T,double *virial);
double sample_binodal(double vmax);
void euler_step(double *x, double *y, double *z,
                double *vx, double *vy, double *vz,
                double *fx, double *fy, double *fz,double L, double dt,double *virial);
void thermostat(double *vx, double *vy, double *vz, double T, double dt);
int main() {
    srand(time(NULL));

    double x[N], y[N], z[N];
    double xo[N], yo[N], zo[N];
    double vx[N], vy[N], vz[N];
    double fx[N], fy[N], fz[N];

    // --- DECLARACIONES ADICIONALES PARA MSD CON PBC ---
    double x_unwrapped[N], y_unwrapped[N], z_unwrapped[N]; // posiciones desenvueltas
    double x_prev[N], y_prev[N], z_prev[N];               // posiciones previas para acumulación

    double dts[] = {0.001};
    int ndt = sizeof(dts)/sizeof(dts[0]);
    int irho; 
    double rho;
    double T = 1.1; // temperatura en unidades reducidas

    double L; // longitud de la caja
    int n_cells = ceil(cbrt(N));

    char filename[64];

    for(irho=0; irho<nrho; irho++) {
        rho = rhos[irho];
        L = cbrt(N / rho);

        for(int idt=0; idt<ndt; idt++) {
            double dt = dts[idt];
            double virial;

            double t_total = 1000.0;
            int steps = (int)(t_total / dt);

            // Inicialización de posiciones y velocidades
            init_sc(x, y, z, n_cells, L);
            for(int i=0; i<N; i++){
                xo[i] = x[i];
                yo[i] = y[i];
                zo[i] = z[i];

                x_unwrapped[i] = x[i];
                y_unwrapped[i] = y[i];
                z_unwrapped[i] = z[i];

                x_prev[i] = x[i];
                y_prev[i] = y[i];
                z_prev[i] = z[i];
            }

            init_velocities(vx, vy, vz, T);
            compute_LJ_forces_pbc(x, y, z, fx, fy, fz, L, &virial);

            sprintf(filename, "msd.dat");
            FILE *f = fopen(filename, "w");
            if(!f){ perror("Error al abrir archivo"); return 1; }
            fprintf(f, "# step time MSD\n");

            for(int step=0; step<steps; step++){
                // Avanzar la simulación
                verlet_velocity_step_thermostat(x, y, z, vx, vy, vz, fx, fy, fz, L, dt, T, &virial);

                // Actualizar posiciones desenvueltas
                for(int i=0; i<N; i++){
                    double dx_step = minimum_image(x[i] - x_prev[i], L);
                    double dy_step = minimum_image(y[i] - y_prev[i], L);
                    double dz_step = minimum_image(z[i] - z_prev[i], L);

                    x_unwrapped[i] += dx_step;
                    y_unwrapped[i] += dy_step;
                    z_unwrapped[i] += dz_step;

                    x_prev[i] = x[i];
                    y_prev[i] = y[i];
                    z_prev[i] = z[i];
                }

                if(step % 10 == 0){
                    // Calcular MSD usando posiciones desenvueltas
                    double msd = 0.0;
                    for(int i=0; i<N; i++){
                        double dx = x_unwrapped[i] - xo[i];
                        double dy = y_unwrapped[i] - yo[i];
                        double dz = z_unwrapped[i] - zo[i];
                        msd += dx*dx + dy*dy + dz*dz;
                    }
                    msd /= N;

                    double KE = kinetic_energy(vx, vy, vz);
                    double PE = compute_LJ_potential(x, y, z, L);
                    double TE = KE + PE;

                    double T_inst = (2.0*KE)/(3.0*N);
                    double V = L*L*L;
                    double P = rho*T_inst + virial/(3.0*V);

                    fprintf(f, "%d %.6f %.6f\n", step, step*dt, msd);
                }
            }
            fclose(f);
            printf("Archivo %s generado.\n", filename);
        }
    }

    return 0;
}



// --- Funciones auxiliares ---

static inline double SafeRandom(){
    double u = (double)rand()/RAND_MAX;
    if(u<=0.0) u=1e-16;
    if(u>=1.0) u=1.0-1e-16;
    return u;
}

void BoxMuller(double* num1, double* num2){
    double u1 = SafeRandom();
    double u2 = SafeRandom();
    double r = sqrt(-2.0*log(u1));
    double theta = 2.0*pi*u2;
    *num1 = r*cos(theta);
    *num2 = r*sin(theta);
}

void init_sc(double *x, double *y, double *z, int n_cells, double L){
    double a = L / n_cells;
    int idx=0;
    for(int i=0;i<n_cells;i++)
        for(int j=0;j<n_cells;j++)
            for(int k=0;k<n_cells;k++)
                if(idx<N){
                    x[idx] = i*a;
                    y[idx] = j*a;
                    z[idx] = k*a;
                    idx++;
                }
}


double sample_binodal(double vmax){
    return (rand() % 2) ? vmax : -vmax;
}

void init_velocities(double *vx, double *vy, double *vz, double T){
    double v_max = sqrt(KB*T/mass);  // desviación estándar para la distribución

    // Generar velocidades aleatorias binodales
    for(int i=0; i<N; i++){
        vx[i] = sample_binodal(v_max);
        vy[i] = sample_binodal(v_max);
        vz[i] = sample_binodal(v_max);
    }

    // Quitar momento total
    double px=0, py=0, pz=0;
    for(int i=0; i<N; i++){
        px += vx[i]; py += vy[i]; pz += vz[i];
    }
    px /= N; py /= N; pz /= N;
    for(int i=0; i<N; i++){
        vx[i] -= px; vy[i] -= py; vz[i] -= pz;
    }
}


void save_xyz(const char *filename, double *x, double *y, double *z){
    FILE *f = fopen(filename,"w");
    if(!f){ perror("Error abriendo XYZ"); return; }
    fprintf(f,"%d\n\n", N);
    for(int i=0;i<N;i++)
        fprintf(f,"Ar %.6f %.6f %.6f\n", x[i], y[i], z[i]);
    fclose(f);
}

void save_velocities(const char *filename, double *vx, double *vy, double *vz){
    FILE *f = fopen(filename,"w");
    if(!f){ perror("Error abriendo velocidades"); return; }

    fprintf(f,"# vx vy vz |v|\n");
    for(int i=0;i<N;i++){
        double vmod = sqrt(vx[i]*vx[i] +vy[i]*vy[i] +vz[i]*vz[i]);
        fprintf(f,"%.6f %.6f %.6f %.6f\n",vx[i], vy[i], vz[i], vmod);
    }
    fclose(f);
}


double minimum_image(double dx, double L){
    if(dx > 0.5*L) dx -= L;
    if(dx < -0.5*L) dx += L;
    return dx;
}

void apply_pbc_positions(double *x, double *y, double *z, double L){
    for(int i=0;i<N;i++){
        if(x[i]>=L) x[i]-=L; if(x[i]<0) x[i]+=L;
        if(y[i]>=L) y[i]-=L; if(y[i]<0) y[i]+=L;
        if(z[i]>=L) z[i]-=L; if(z[i]<0) z[i]+=L;
    }
}

void compute_LJ_forces_pbc(double *x, double *y, double *z,
                           double *fx, double *fy, double *fz,
                            double L,double *virial){

    double cutoff2 = cutoff*cutoff;
    const double r2_min = 1e-12;

    for(int i=0;i<N;i++){
        fx[i]=fy[i]=fz[i]=0.0;
    }
    (*virial) = 0.0; 

    for(int i=0;i<N-1;i++){
        for(int j=i+1;j<N;j++){
            double dx = minimum_image(x[i]-x[j], L);
            double dy = minimum_image(y[i]-y[j], L);
            double dz = minimum_image(z[i]-z[j], L);

            double r2 = dx*dx + dy*dy + dz*dz;

            if(r2<cutoff2 && r2>r2_min){
                double r2i = 1.0/r2;
                double r6i = r2i*r2i*r2i;
                double f = 48.0*r6i*(r6i-0.5)*r2i;

                double fxij = f*dx;
                double fyij = f*dy;
                double fzij = f*dz;

                fx[i]+=fxij; fy[i]+=fyij; fz[i]+=fzij;
                fx[j]-=fxij; fy[j]-=fyij; fz[j]-=fzij;

                // virial
                (*virial) += dx*fxij + dy*fyij + dz*fzij;
            }
        }
    }
}

double compute_LJ_potential(double *x, double *y, double *z,
                             double L){
    double E=0.0;
    double cutoff2=cutoff*cutoff;
    for(int i=0;i<N-1;i++)
        for(int j=i+1;j<N;j++){
            double dx = minimum_image(x[i]-x[j], L);
            double dy = minimum_image(y[i]-y[j], L);
            double dz = minimum_image(z[i]-z[j], L);
            double r2 = dx*dx + dy*dy + dz*dz;
            if(r2<cutoff2){
                double r2i=1.0/r2;
                double r6i=r2i*r2i*r2i;
                E += 4.0*r6i*(r6i-1.0);
            }
        }
    return E;
}

void euler_step( double *x, double *y, double *z,
                double *vx, double *vy, double *vz,
                double *fx, double *fy, double *fz,
              double L, double dt,double *virial)
{
    for (int i = 0; i < N; i++) {
        vx[i] += fx[i] * dt;
        vy[i] += fy[i] * dt;
        vz[i] += fz[i] * dt;

        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
        z[i] += vz[i] * dt;
    }

    apply_pbc_positions(x, y, z, L);
    compute_LJ_forces_pbc(x, y, z, fx, fy, fz,L,virial);
}

double kinetic_energy(double *vx, double *vy, double *vz){
    double E=0.0;
    for(int i=0;i<N;i++) E+=0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    return E;
}

void verlet_velocity_step(double *x, double *y, double *z,
                          double *vx, double *vy, double *vz,
                          double *fx, double *fy, double *fz,
                          double L, double dt,double *virial){
    for(int i=0;i<N;i++){
        vx[i]+=0.5*fx[i]*dt;
        vy[i]+=0.5*fy[i]*dt;
        vz[i]+=0.5*fz[i]*dt;
        x[i]+=vx[i]*dt;
        y[i]+=vy[i]*dt;
        z[i]+=vz[i]*dt;
    }

    apply_pbc_positions(x,y,z,L);
    compute_LJ_forces_pbc(x,y,z,fx,fy,fz,L,virial);

    for(int i=0;i<N;i++){
        vx[i]+=0.5*fx[i]*dt;
        vy[i]+=0.5*fy[i]*dt;
        vz[i]+=0.5*fz[i]*dt;
    }
}

void verlet_velocity_step_thermostat(double *x, double *y, double *z,
                          double *vx, double *vy, double *vz,
                          double *fx, double *fy, double *fz,
                          double L, double dt, double T,double *virial){
    for(int i=0;i<N;i++){
        vx[i]+=0.5*fx[i]*dt;
        vy[i]+=0.5*fy[i]*dt;
        vz[i]+=0.5*fz[i]*dt;
        x[i]+=vx[i]*dt;
        y[i]+=vy[i]*dt;
        z[i]+=vz[i]*dt;
    }

    apply_pbc_positions(x,y,z,L);
    compute_LJ_forces_pbc(x,y,z,fx,fy,fz,L,virial);

    for(int i=0;i<N;i++){
        vx[i]+=0.5*fx[i]*dt;
        vy[i]+=0.5*fy[i]*dt;
        vz[i]+=0.5*fz[i]*dt;
    }

    thermostat(vx,vy,vz,T,dt);
}

void thermostat(double *vx, double *vy, double *vz, double T, double dt){
    double p = nu * dt;

    double z1,z2,z3,z4;

    for (int i = 0; i < N; i++) {
        BoxMuller(&z1,&z2);
        BoxMuller(&z3,&z4);

        if (SafeRandom() < p) {
            double sigma = sqrt(KB*T/mass);
            vx[i] = sigma*z1;
            vy[i] = sigma*z2;
            vz[i] = sigma*z3;
        }
    }
}


