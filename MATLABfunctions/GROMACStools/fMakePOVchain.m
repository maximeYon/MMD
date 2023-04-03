function res = fMakePOVchain(fnam,gro,molnr,chainatoms)

fid = fopen([fnam '.inc'],'w');

for ncount = 1:length(molnr)
    n = molnr(ncount);
    index = [];
    for m = 1:length(chainatoms)
        index = [index; find(all([strcmp(gro.atom, chainatoms{m}) gro.molnr==n],2))];
    end
    
%     thresh = .5;
%     for m = 2:length(index)
% %             b4 = (gro.z(index(m))-gro.z(index(m-1)));
% %             moveflags = [0 0];
%         if (gro.x(index(m))-gro.x(index(m-1))) > thresh*gro.boxxx
%             gro.x(index(m)) = gro.x(index(m)) - gro.boxxx;
%         end
%         if (gro.x(index(m))-gro.x(index(m-1))) < -thresh*gro.boxxx
%             gro.x(index(m)) = gro.x(index(m)) + gro.boxxx;
%         end
%         if (gro.y(index(m))-gro.y(index(m-1))) > thresh*gro.boxyy
%             gro.y(index(m)) = gro.y(index(m)) - gro.boxyy;
%         end
%         if (gro.y(index(m))-gro.y(index(m-1))) < -thresh*gro.boxyy
%             gro.y(index(m)) = gro.y(index(m)) + gro.boxyy;
%         end
%         if (gro.z(index(m))-gro.z(index(m-1))) > thresh*gro.boxzz
%             gro.z(index(m)) = gro.z(index(m)) - gro.boxzz;
%             %moveflags(1) = 1;
%         end
%         if (gro.z(index(m))-gro.z(index(m-1))) < -thresh*gro.boxzz
%             gro.z(index(m)) = gro.z(index(m)) + gro.boxzz;
%             %moveflags(2) = 1;
%         end
% %             af = (gro.z(index(m))-gro.z(index(m-1)));
% %             [b4 af moveflags]
% %     figure(1), clf
% %     plot3(gro.x(index),gro.y(index),gro.z(index),'-',gro.x(index(m)),gro.y(index(m)),gro.z(index(m)),'o')
% %     axis(.75*gro.boxzz*[-1 1 -1 1 -1 1])
% %     axis square    
% %     view(0,0)
% %     title([num2str(n) '   ' num2str(m)])
% %     pause%(.5)
%     end
    
% %     indexmove = find((gro.x(index)-gro.x(index(indexmid))) < -thresh*gro.boxxx);
%     gro.x(index(indexmove)) = gro.x(index(indexmove)) + gro.boxxx;
%     indexmove = find((gro.y(index)-gro.y(index(indexmid))) > thresh*gro.boxyy);
%     gro.y(index(indexmove)) = gro.y(index(indexmove)) - gro.boxyy;
%     indexmove = find((gro.y(index)-gro.y(index(indexmid))) < -thresh*gro.boxyy);
%     gro.y(index(indexmove)) = gro.y(index(indexmove)) + gro.boxyy;
%     indexmove = find((gro.z(index)-gro.z(index(indexmid))) > thresh*gro.boxzz);
%     gro.z(index(indexmove)) = gro.z(index(indexmove)) - gro.boxzz;
%     indexmove = find((gro.z(index)-gro.z(index(indexmid))) < -thresh*gro.boxzz);
%     gro.z(index(indexmove)) = gro.z(index(indexmove)) + gro.boxzz;
%     indexmove = find((gro.y(index)-min(gro.y(index))) > .5*gro.boxyy);
%     gro.y(index(indexmove)) = gro.y(index(indexmove)) - gro.boxyy;
%     indexmove = find((gro.z(index)-min(gro.z(index))) > .5*gro.boxzz);
%     gro.z(index(indexmove)) = gro.z(index(indexmove)) - gro.boxzz;
    
    tempstr = ['sphere{<' num2str(gro.x(index(1))) ','...
        num2str(gro.z(index(1))) ',' ...
        num2str(gro.y(index(1))) '>' ',rad}']; 
    fprintf(fid,'%s\n',tempstr);
    for m = 2:length(index)
        tempstr = ['cylinder{<' num2str(gro.x(index(m-1))) ','...
            num2str(gro.z(index(m-1))) ',' ...
            num2str(gro.y(index(m-1))) '>,<'...
            num2str(gro.x(index(m))) ','...
            num2str(gro.z(index(m))) ',' ...
            num2str(gro.y(index(m))) '>,rad}']; 
        fprintf(fid,'%s\n',tempstr);
        tempstr = ['sphere{<' num2str(gro.x(index(m))) ','...
            num2str(gro.z(index(m))) ',' ...
            num2str(gro.y(index(m))) '>' ',rad}']; 
        fprintf(fid,'%s\n',tempstr);
    end
end
fclose(fid);

res = 1;
