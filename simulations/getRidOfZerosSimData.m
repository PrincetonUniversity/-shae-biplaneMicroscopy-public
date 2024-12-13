function traj = getRidOfZerosSimData (traj)

temp = traj(:,1);
for iiL= 1:length(traj(:,1))
    if isnan(traj(iiL, 1))
        if iiL ==1
            before = nan;
        else
            before = traj(iiL-1,1);
        end
        if iiL == length (traj(:,1))
            after= nan;
        else
            after = traj(iiL+1,1);
        end
        
        
        cMin =0;
        cMax =0;
        while isnan(before)
            before = traj(iiL-cMin,1);
            cMin =cMin+1;
            if (iiL-cMin==0)
                break
            end
        end
        while isnan(after)
            after = traj(iiL+cMax,1);
            cMax = cMax +1;
            if (iiL+cMax>=length(traj))
                break
            end
        end
        if iiL==1
            temp(iiL)=after;
        elseif iiL ==length(traj(:,1))
            temp(iiL)=before;
        else 
            temp(iiL)=(before+after)/2;
        end
    end
end
traj(:,1)=temp;

end