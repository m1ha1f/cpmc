function [F,D] = csift(I, spacing, scales, color_descriptor_type)
    if(size(I,3) == 3)
        twist = floor(sum(sum(sum(I))));
    else
        twist = floor(sum(sum(I)));
    end

    DefaultVal('*color_descriptor_type', '''csift''');

    rand('twister', twist);

    img_filename = ['./img' int2str(rand*10000000000000) '.jpg'];
    imwrite(I, img_filename);

    if(size(I,3) == 3)
        twist = floor(sum(sum(sum(I))));
    else
        twist = floor(sum(sum(I)));
    end

    rand('twister', twist);

    output_filename = ['./csift' int2str(rand*10000000000000) '.tmp'];
    path = which('csift');
    path = path(1:end-length('csift.m'));

    scales_str = [];
    for i=1:length(scales)
        scales_str = [scales_str sprintf('%f+', scales(i))];
    end
    scales_str(end) = [];

    if strcmp(computer,'GLNX86')
        exec_name='colorDescriptor32';
    else
        exec_name='colorDescriptor64';
    end
    system([path exec_name ' ' img_filename  ' --detector densesampling --ds_spacing ' int2str(spacing) ' --ds_scales ' scales_str ' --descriptor ' color_descriptor_type ' --output ' output_filename]);

    %system([path 'colorDescriptor64 '  img_filename  ' --detector densesampling --ds_spacing ' int2str(spacing) ' --ds_scales ' scales_str ' --descriptor ' color_descriptor_type ' --output ' output_filename]);
    %system(['grep CIRCLE ' output_filename ' > ' output_filename '2']);
    fid = fopen(output_filename);
    [n_feats] = textscan(fid, 'KOEN1\n%d\n%u');
    n_feats = n_feats{1};
    [data] = textscan(fid, '%s', 'delimiter', ';');
    nrows = size(data{1},1);
    D = cell(n_feats,1);
    F = cell(n_feats, 1);
    for i=2:2:nrows
        D{i/2} = str2num(data{1}{i});
    end
    
    D = cell2mat(D)';
    if(strcmp(color_descriptor_type, 'rgbhistogram') || strcmp(color_descriptor_type, 'huehistogram'))        
       D = single(D);
    else
        D = uint8(D);
    end
        
    
    for i=1:2:nrows
        F{floor((i/2))+1} = str2num(data{1}{i}(8:end-1));
    end
    F = cell2mat(F)';
    fclose(fid);
    delete(output_filename);
    %delete([output_filename '2']);
    delete(img_filename);
    F = single(F);
end