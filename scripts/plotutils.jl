function add_panel_letters!(panels; fontsize=9)
    letters = 'A':'Z'
    for (p, letter) in zip(panels, letters)
        annotate!(p, (0.05, 0.95), text("($letter)", :left, fontsize))
    end
end