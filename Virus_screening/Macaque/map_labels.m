function new_label = get_new_label(old_label, original_labels, new_labels)
    idx = strcmp(old_label, original_labels);
    if any(idx)
        new_label = new_labels{idx};
    else
        new_label = old_label;  % Keep original if no match found
    end
end