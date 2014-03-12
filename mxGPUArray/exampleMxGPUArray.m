

x = gpuArray.ones(4,4);
y = mexGPUExample(x)


disp(['class(x) = ',class(x),', class(y) = ',class(y)])
