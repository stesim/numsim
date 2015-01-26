function animate(ranks,instance)
    filename = ['field-uvpT-instance-' num2str(instance,'%03d') '-step-%d-rank-000-000.txt#'];
    list = dir(fullfile(cd, 'field-uvpT-instance-*-step-*-rank-*-*.txt'));
    name = {list.name};
    str  = sprintf('%s#', name{:});
    num  = sscanf(str, filename);
    [dummy, index] = sort(num);
    %name = name(index);
    
    for i=1:length(num)
        visualize(ranks, num(i), instance);
        pause(0.1);
    end
end