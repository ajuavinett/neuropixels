function clearObjects

availObjects = instrfind;
if length(availObjects) > 1;
    for iObject = 2:length(availObjects);
        extraObject = availObjects(iObject);
        fclose(extraObject);
        delete(extraObject);
    end
end

disp(strcat('SuccessfullyCleared ',num2str(length(availObjects)), 'object(s)'))
