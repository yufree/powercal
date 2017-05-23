function [SIGMAlow]=SurfacePlot(output, variable,metric,correction, sizeeff,samplsizes,nreps)


MUtot = output{1,variable}{correction,metric};
[NS, NSE]=size(MUtot);
SIGMAtot = output{1,variable}{correction,metric+5};

SIGMAlow=MUtot-1.96*SIGMAtot/sqrt(nreps);

for i=1:NS
    for j=1:NSE
        if(SIGMAlow(i,j)<0)
            SIGMAlow(i,j)=0;
        end
    end
end

%plot options%

grey=gray;
grey(20:24,:)=[];
grey(1:4,:)=[];

minSIGMAlow=min(min(SIGMAlow));
scaledimg = (floor(((SIGMAlow - minSIGMAlow) ./ (max(max(SIGMAlow)) - minSIGMAlow)) * 255));
colorimg = ind2rgb(scaledimg,jet(256));

surf(sizeeff,samplsizes,MUtot,'edgecolor','none');
colormap(grey)
alpha(0.9)
hold on

surf([min(sizeeff) max(sizeeff)],[min(samplsizes) max(samplsizes)],repmat(-0.6, [2 2]),colorimg,'facecolor','texture')

xlabel ('Sample size','Fontsize',12);
ylabel ('Effect Size','Fontsize',12);
zlabel ('Rate','Fontsize',12);
