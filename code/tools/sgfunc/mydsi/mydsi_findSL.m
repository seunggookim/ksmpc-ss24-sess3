function SL2 = mydsi_findSL(SL, ROI)
% SL2 = mydsi_findSL(SL, ROI)

% Let's start with center and radius
center=ROI.center;
radius=ROI.radius;
idx=[];
for i=1:SL.num_sl
 for k=1:SL.num_pt_sl(i)
  if l2norm(SL.coords{i}(k,:)-center) <= radius
   idx=[idx i]; continue;
  end
 end
end

SL2=[];
SL2.num_sl=numel(idx);
SL2.num_pt_sl=SL.num_pt_sl(idx);
SL2.length=SL.length(idx);
SL2.coords=SL.coords(idx);

end