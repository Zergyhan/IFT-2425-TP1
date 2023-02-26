//
// Felix Rouleau; 20188552; felix.rouleau@gmail.com
//
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <wait.h>
#include <malloc.h>
#include <stdbool.h>

enum derivativeTypes {
    QUESTION_A,
    QUESTION_B,
    QUESTION_C
};

double newtonApproximation(double cmv, double *y, int lenY) {
    // Function is sum from 1->N (yi^cmv*ln(yi)) / sum 1->N yi^cmv - 1/cmv * "constant"
    // "Constant" is 1/N sum 1->N ln(yi).
    double total_num = 0;
    double total_denum = 0;
    double rightHand = 0;
    double total = 0;
    for (int i = 0; i > lenY; i++) {
        total_num += pow(y[i], cmv) * log(y[i]);
        total_denum += pow(y[i], cmv);
        rightHand += log(y[i]);
    }
    total = total_num / total_denum;
    total -= 1 / cmv;
    total -= rightHand / lenY;
    return total;
}

double newtonDerivative(double cmv, double *y, int lenY, enum derivativeTypes type) {
    switch (type) {
        case QUESTION_A:
            // Version where it is f'(x) = [f(x+h) - f(x)] / h where h is 10e-6
            return (newtonApproximation(cmv + 10e-6, y, lenY) - newtonApproximation(cmv, y, lenY)) / 10e-6;
        case QUESTION_B:
            // Version where it is f'(x) = [-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)] / 12h where h is 10e-6
            return (-newtonApproximation(cmv + 2 * 10e-6, y, lenY) +
                    8 * newtonApproximation(cmv + 10e-6, y, lenY) -
                    8 * newtonApproximation(cmv - 10e-6, y, lenY) +
                    newtonApproximation(cmv - 2 * 10e-6, y, lenY)) / (12 * 10e-6);
        case QUESTION_C:
            // Version where it is the analytic derivative of f(x)
            break;
    }
}

double cmvEstimator(double *y, int lenY, enum derivativeTypes type) {
    double tolerance = 10e-6;
    double cmv = 0.25;
    double old_cmv;
    do {
        old_cmv = cmv;
        cmv = cmv - newtonApproximation(cmv, y, lenY) / newtonDerivative(cmv, y, lenY, type);
    } while (fabs(cmv - old_cmv) > tolerance &&
             newtonApproximation(cmv, y, lenY) > tolerance &&
             newtonDerivative(cmv, y, lenY, type) != 0);
}

int WIDTH = 512;
int HEIGHT = 512;

typedef struct ComplexPoint {
    long double x;
    long double y;
} ComplexPoint;

/// \brief Add two complex numbers
/// \param a ComplexPoint a to add
/// \param b ComplexPoint b to add
/// \return ComplexPoint a + b
ComplexPoint addComplex(ComplexPoint a, ComplexPoint b) {
    ComplexPoint result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    return result;
}

/// \brief Square a complex number
/// \param z Complex number to square
/// \return ComplexPoint z^2
ComplexPoint squareComplex(ComplexPoint z) {
    ComplexPoint result;
    result.x = z.x * z.x - z.y * z.y;
    result.y = 2 * z.x * z.y;
    return result;
}

long double iterateMandelbrot(long double x, long double y, int maxIterations) {
    ComplexPoint z = {0, 0};
    ComplexPoint c = {x, y};
    for (int iteration = 0; iteration < maxIterations; iteration++) {
        z = addComplex(squareComplex(z), c);
        // If the magnitude of z is greater than 2, then it is not in the Mandelbrot set
        if (z.x * z.x + z.y * z.y > 4) {
            // Return the number from 0 to 1 that represents how close it is to the Mandelbrot set
            return iteration / (long double) maxIterations;
        }
    }
    return 1.0;
}

/// \brief Get the pixel from a point in the complex plane
/// \param x Point in the complex plane
/// \param y Point in the complex plane
/// \param pixel[2] Output array of size 2, [0] is x, [1] is y
void getPixelFromPoint(long double x, long double y, int pixel[2]) {
    pixel[0] = (int) floorl(x * (WIDTH - 1) / 2 + WIDTH / 1.35);
    pixel[1] = (int) floorl(y * (HEIGHT - 1) / 2.0 + HEIGHT / 2.0);
}

/// \brief Get the history of a point in the complex plane if
/// \param horizontal
/// \param vertical
/// \param maxIterations
/// \return NULL if the point is in the Mandelbrot set, otherwise a pointer to a history of ComplexPoints
ComplexPoint *
iterateHiddenMandelbrot(long double horizontal, long double vertical, int maxIterations, bool isInterior) {
    ComplexPoint *history = malloc(sizeof(ComplexPoint) * maxIterations);
    ComplexPoint z = {0, 0};
    ComplexPoint c = {horizontal, vertical};
    bool escapes = false;

    for (int iteration = 0; iteration < maxIterations; iteration++) {
        z = addComplex(squareComplex(z), c);
        history[iteration] = z;
        // If the magnitude of z is greater than 2, then it is not in the Mandelbrot set
        if (z.x * z.x + z.y * z.y > 4) {
            escapes = true;
        }
    }
    if (isInterior) {
        if (escapes) return (ComplexPoint *) NULL;
        return history;

    } else {
        if (escapes) return history;
        return (ComplexPoint *) NULL;
    }
}

/// \brief Create a Mandelbrot set image
/// \param maxIterations Maximum number of iterations to calculate z = z^2 + c
void createMandelbrot(int maxIterations) {
    // Question 2.1
    long double image[HEIGHT][WIDTH];
    for (int i = 0; i < HEIGHT; i++) {
        for (int j = 0; j < WIDTH; j++) {
            image[i][j] = 0;
        }
    }
    FILE *fp = fopen("image.ppm", "wb");
    fprintf(fp, "P6\n%d %d\n255\n", HEIGHT, WIDTH);

    clock_t start = clock();

    for (int l = 0; l < HEIGHT; l++) {
        for (int k = 0; k < WIDTH; k++) {
            long double colour = iterateMandelbrot(2 * (k - WIDTH / 1.35) / (WIDTH - 1),
                                                   2 * (l - HEIGHT / 2.0) / (HEIGHT - 1),
                                                   maxIterations);
            int finalColour = 255 * (1 - (int) floorl(colour));
            unsigned char writeColour[3] = {finalColour, finalColour, finalColour};
            fwrite(writeColour, 1, 3, fp);
        }
    }

    fclose(fp);
    clock_t total = clock() - start;
    printf("Time taken for question 2.1: %f\n", (double) total / CLOCKS_PER_SEC);

    // Fork and open image with xdg-open
    pid_t pid = fork();
    if (pid == 0) {
        execlp("xdg-open", "xdg-open", "image.ppm", NULL);
    } else {
        // Wait for the image viewer to close
        wait(NULL);
    }
}

/// \brief Create a hidden Mandelbrot set image
/// \param maxIterations Maximum number of iterations to calculate z = z^2 + c
void createHiddenMandelbrot(int maxIterations, bool isInterior, int PIXEL_INTENSITY) {
    // Question 2.2
    // Take a point k, l. If the point c escapes, then add +1 to the pixels where the z travelled to.
    long double image[HEIGHT][WIDTH];
    for (int i = 0; i < HEIGHT; i++) {
        for (int j = 0; j < WIDTH; j++) {
            image[i][j] = 0;
        }
    }

    clock_t start = clock();
    clock_t last_iteration = start;
    // Step is 1/1000th
    int STEP = 2;
    int numToIterate = 1000;
    for (int vertical = 0; vertical <= numToIterate; vertical++) {
        for (int horizontal = 0; horizontal <= numToIterate; horizontal++) {
            long double x = (long double) (horizontal * STEP) / 1000 - 1.5;
            long double y = (long double) (vertical * STEP) / 1000 - 1.0;
            ComplexPoint *history = iterateHiddenMandelbrot(x, y, maxIterations, isInterior);
            if (history == NULL) continue;
            for (int i = 0; i < maxIterations; i++) {
                int pixel[2];
                getPixelFromPoint(history[i].x, history[i].y, pixel);
                if (pixel[0] < 0 || pixel[1] < 0
                    || pixel[0] >= WIDTH || pixel[1] >= HEIGHT)
                    continue;

                image[pixel[1]][pixel[0]] += PIXEL_INTENSITY;
            }
            free(history);
        }
        if (vertical % 50 == 0) {
            if (vertical == 0) continue;
            printf("Iteration %d/%d\nTime taken on average:%lds\n", vertical, numToIterate,
                   ((clock() - last_iteration) / CLOCKS_PER_SEC));
            printf("Estimated time remaining: %lds\n\n",
                   ((clock() - last_iteration) / CLOCKS_PER_SEC) * (numToIterate - vertical) / 50);
            last_iteration = clock();
        }
    }

    char *filename;
    if (isInterior) {
        filename = "image3.ppm";
    } else {
        filename = "image2.ppm";
    }
    FILE *fp = fopen(filename, "wb");
    fprintf(fp, "P6\n%d %d\n255\n", HEIGHT, WIDTH);
    for (int l = 0; l < HEIGHT; l++) {
        for (int k = 0; k < WIDTH; k++) {
            int finalColour = (int) image[l][k];
            if (finalColour > 255) finalColour = 255;
            unsigned char writeColour[3] = {finalColour, finalColour, finalColour};
            fwrite(writeColour, 1, 3, fp);
        }
    }
    fclose(fp);
    clock_t total = clock() - start;
    printf("Time taken for question 2.2: %f\n", (double) total / CLOCKS_PER_SEC);

    // Fork and open image with xdg-open
    pid_t pid = fork();
    if (pid == 0) {
        execlp("xdg-open", "xdg-open", filename, NULL);
    } else {
        // Wait for the image viewer to close
        wait(NULL);
    }
}

int main() {
    int question = 1;
    if (question == 0) {
        double y[] = {0.11, 0.24, 0.27, 0.52, 1.13, 1.54, 1.71, 1.84, 1.92, 2.01};
        int lenY = 10;
        double cmv = cmvEstimator(y, lenY, QUESTION_A);
        printf("CMV is %f", cmv);
    } else if (question == 1) {
        int maxIterations = 200;
        createMandelbrot(maxIterations);
        createHiddenMandelbrot(maxIterations, false, 5);
        createHiddenMandelbrot(maxIterations, true, 1);

    }
    return 0;
}
