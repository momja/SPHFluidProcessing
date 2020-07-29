boolean[][] binarizeImage(String file, int threshMax) {
    PImage img = loadImage(file);

    float rWeight = 0.2989;
    float gWeight = 0.5870;
    float bWeight = 0.1140;

    boolean[][] bitmap = new boolean[img.height][img.width];

    img.loadPixels();
    for (int i = 0; i < img.height; i++) {
        for (int j = 0; j < img.width; j++) {
            color pix = img.pixels[j+img.width*i];
            if (rWeight*red(pix) + gWeight*green(pix) + bWeight*blue(pix) < threshMax) {
                bitmap[i][j] = true;
            } else {
                bitmap[i][j] = false;
            }
        }
    }
    return bitmap;
}