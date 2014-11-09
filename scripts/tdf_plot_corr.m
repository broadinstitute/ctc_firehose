%%%%%%%%%%%%%%%%%%%%
%Single cell QC:
%Generate image that combines the autocorrelation coefficient and image across all samples in set
%
%MR 9 June 2014
%%%%%%%%%%%%%%%%%%%

%using matlab:

function tdf_plot_corr(sample, image_path, corr)

	display(['plotting for sample: ', sample])
	I = imread(image_path);

	plot_text = ['corr coefficient: ', corr];
	image(I)
	text(400,300,plot_text)
	text(400,350, sample)
	axis off
	print_D(sample,{{'pdf'},{'png','-r180'}});
end