% this script is to experiment reading in the brain images

imDir = 'C:\Users\tnair_000\Documents\medimage\OASIS\disc1\OAS1_0001_MR1\PROCESSED\MPRAGE\T88_111';
hdr = [imDir '\OAS1_0001_MR1_mpr_n4_anon_111_t88_masked_gfc.hdr'];

hdrInfo = analyze75info(hdr);
X = analyze75read( hdrInfo );
slice = im2double(X(:,:,100));
imagesc(slice);