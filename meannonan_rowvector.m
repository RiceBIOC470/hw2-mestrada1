function mm=meannonan_rowvector(x)
for ii = 1:length(x)
    for y = 1:size(x,1)
        if isnan(x(y, ii)) == 1
        x(y, ii) = 0;
        end
    end
end
denom = zeros(1,length(x));
for ii = 1:length(x)
    dm = size(x,1);
    for y = 1:size(x,1)
        if x(y,ii) == 0
        dm = dm - 1;
        end
    end 
    denom(ii) = dm;
end
meannonan = zeros(1,5);
for ii = 1:length(x)
    meannonan(ii) = sum(x(:,ii))/(denom(ii));
end
mm=meannonan;