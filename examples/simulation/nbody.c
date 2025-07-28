#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Constants
const double G = 6.67430e-11; // Gravitational constant
const double DT = 0.01;       // Time step

// Body struct
typedef struct {
    double x, y;       // Position
    double vx, vy;     // Velocity
    double ax, ay;     // Acceleration
    double mass;
} Body;

// Function to compute gravitational force
void compute_forces(Body* bodies, int n) {
    for (int i = 0; i < n; i++) {
        bodies[i].ax = 0;
        bodies[i].ay = 0;
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double dist_sq = dx*dx + dy*dy + 1e-10; // Softening
            double dist = sqrt(dist_sq);
            double force = G * bodies[j].mass / dist_sq;
            bodies[i].ax += force * dx / dist;
            bodies[i].ay += force * dy / dist;
        }
    }
}

// Function to update positions using simple Euler integration (for simplicity here)
void update_bodies(Body* bodies, int n) {
    for (int i = 0; i < n; i++) {
        bodies[i].vx += bodies[i].ax * DT;
        bodies[i].vy += bodies[i].ay * DT;
        bodies[i].x += bodies[i].vx * DT;
        bodies[i].y += bodies[i].vy * DT;
    }
}

// Function to initialize bodies randomly
void init_bodies(Body* bodies, int n) {
    for (int i = 0; i < n; i++) {
        bodies[i].x = rand() % 1000;
        bodies[i].y = rand() % 1000;
        bodies[i].vx = 0;
        bodies[i].vy = 0;
        bodies[i].mass = 1e5 + rand() % 10000;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        printf("Usage: %s <num_bodies> <num_steps>\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);       // Number of bodies
    int steps = atoi(argv[2]);   // Number of time steps

    Body* bodies = (Body*)malloc(n * sizeof(Body));
    init_bodies(bodies, n);

    clock_t start = clock();

    for (int step = 0; step < steps; step++) {
        compute_forces(bodies, n);
        update_bodies(bodies, n);
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Simulation completed in %.4f seconds\n", elapsed);

    // Optional: print final state of first few bodies
    for (int i = 0; i < (n < 5 ? n : 5); i++) {
        printf("Body %d: pos(%.2f, %.2f) vel(%.2f, %.2f)\n",
               i, bodies[i].x, bodies[i].y, bodies[i].vx, bodies[i].vy);
    }

    free(bodies);
    return 0;
}
