int kernelSize = 5; 
double sigma = 1.4;
double kernel[5];
double sumKernel = 0.0;

// generate 1d kernel
for (int i = -kernelSize / 2; i <= kernelSize / 2; i++) {
    kernel[i + kernelSize / 2] = exp(-(i * i) / (2 * sigma * sigma)) / (sqrt(2 * M_PI) * sigma);
    sumKernel += kernel[i + kernelSize / 2];
}
// normalize kernel
for (int i = 0; i < kernelSize; i++) {
    kernel[i] /= sumKernel;
}

// apply 1d gaussian on buffer
void apply1DGaussianBlurLinear(JSAMPLE* src, JSAMPLE* dst, int width, int height, int channels, const double* kernel, int kernelSize, int direction) {
    int halfKernel = kernelSize / 2;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            double sum[3] = { 0.0, 0.0, 0.0 }; // assuming rgb
            for (int k = -halfKernel; k <= halfKernel; k++) {
                int sampleX = x, sampleY = y;
                if (direction == 0) { // horizontal blur
                    sampleX = x + k;
                    if (sampleX < 0) sampleX = 0;
                    if (sampleX >= width) sampleX = width - 1;
                }
                else { // vertical blur
                    sampleY = y + k;
                    if (sampleY < 0) sampleY = 0;
                    if (sampleY >= height) sampleY = height - 1;
                }
                JSAMPLE* pixel = src + (sampleY * width + sampleX) * channels;
                for (int c = 0; c < channels; c++) {
                    sum[c] += pixel[c] * kernel[k + halfKernel];
                }
            }
            JSAMPLE* outPixel = dst + (y * width + x)* channels;
            for (int c = 0; c < channels; c++) {
                outPixel[c] = (JSAMPLE)(sum[c]);
            }
        }
    }
}