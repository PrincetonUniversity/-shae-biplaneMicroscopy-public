function traj = getRidOfNans (traj)

temp = traj(:,3);
for iiL= 1:length(traj(:,3))
    if isnan(traj(iiL, 3))
        if iiL ==1
            before = nan;
        else
            before = traj(iiL-1,3);
        end
        if iiL == length (traj(:,3))
            after= nan;
        else
            after = traj(iiL+1,3);
        end
        
        
        cMin =0;
        cMax =0;
        while isnan(before)
            before = traj(iiL-cMin,3);
            cMin =cMin+1;
            if (iiL-cMin==0)
                break
            end
        end
        while isnan(after)
            after = traj(iiL+cMax,3);
            cMax = cMax +1;
            if (iiL+cMax>=length(traj))
                break
            end
        end
        if iiL==1
            temp(iiL)=after;
        elseif iiL ==length(traj(:,3))
            temp(iiL)=before;
        else 
            temp(iiL)=(before+after)/2;
        end
    end
end
traj(:,3)=temp;

end