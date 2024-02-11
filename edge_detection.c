// libjpeg functions to read jpeg
// convert to grayscale using luminance
// return image data
// parameters:
//   filename: filepath to input file
//   image: pointer to the image data
//   width, height: image dimensions
unsigned char* read_jpeg(char* filename, int* width, int* height) {

    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE* infile;
    JSAMPARRAY pJpegBuffer;
    unsigned char* rowptr;
    unsigned char* raw_image = NULL;

    if ((infile = fopen(filename, "rb")) == NULL) {
        fprintf(stderr, "can't open %s\n", filename);
        return NULL;
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    (void)jpeg_read_header(&cinfo, TRUE);

    (void)jpeg_start_decompress(&cinfo);
    *width = cinfo.output_width;
    *height = cinfo.output_height;

    raw_image = (unsigned char*)malloc(cinfo.output_width * cinfo.output_height);
    rowptr = (unsigned char*)malloc(cinfo.output_width * cinfo.num_components);

    while (cinfo.output_scanline < cinfo.output_height) {
        pJpegBuffer = &rowptr;
        jpeg_read_scanlines(&cinfo, pJpegBuffer, 1);

        for (unsigned int x = 0; x < cinfo.output_width; x++) {
            // convert to grayscale using luminosity
            // TODO: check rowptr for NULL
            int grayVal = 0.299 * rowptr[x * cinfo.num_components] +
                0.587 * rowptr[x * cinfo.num_components + 1] +
                0.114 * rowptr[x * cinfo.num_components + 2];
            // TODO: check raw_image for NULL
            raw_image[cinfo.output_scanline - 1 + x * cinfo.output_width] = (unsigned char)grayVal;
        }
    }

    free(rowptr);
    (void)jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);

    return raw_image;

}

// apply sobel operator to detect edges
// parameters:
//   gray_image: pointer to the image data
//   width, height: image dimensions
//   threshold: above this is an edge
unsigned char* apply_sobel_operator(const unsigned char* gray_image, int width, int height) {
    unsigned char* output_image = (unsigned char*)malloc(width * height);
    if (!output_image) {
        fprintf(stderr, "memory allocation failed for output image\n");
        return NULL;
    }

    int Gx, Gy, sum;

    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            Gx = -gray_image[(y - 1) * width + (x - 1)] - 2 * gray_image[y * width + (x - 1)] - gray_image[(y + 1) * width + (x - 1)]
                + gray_image[(y - 1) * width + (x + 1)] + 2 * gray_image[y * width + (x + 1)] + gray_image[(y + 1) * width + (x + 1)];
            Gy = -gray_image[(y - 1) * width + (x - 1)] - 2 * gray_image[(y - 1) * width + x] - gray_image[(y - 1) * width + (x + 1)]
                + gray_image[(y + 1) * width + (x - 1)] + 2 * gray_image[(y + 1) * width + x] + gray_image[(y + 1) * width + (x + 1)];

            sum = (int)sqrt(Gx * Gx + Gy * Gy);
            output_image[y * width + x] = (sum > 255) ? 255 : sum; // ensure the value in range {0, 255}
        }
    }

    return output_image;
}

// flatten grayscale image pixels to either 0 or 255
// parameters:
//   image: pointer to the image data
//   width, height: image dimensions
//   threshold: above this is an edge
void binarize_image(unsigned char* image, int width, int height, int threshold) {
    for (int i = 0; i < width * height; i++) {
        image[i] = (image[i] > threshold) ? 0 : 255; // edge: black | not edge: white
    }
}
