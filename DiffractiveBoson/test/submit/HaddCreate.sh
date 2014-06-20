mkdir WtoENu_minus
mkdir WtoENu_plus
mkdir WtoMuNu_minus
mkdir WtoMuNu_plus
mkdir ZtoEE_minus
mkdir ZtoEE_plus
mkdir ZtoMuMu_plus
mkdir ZtoMuMu_minus

hadd WtoENu_minus/WtoENu_minus.root PomwigWtoENu_minus/res/*.root
hadd WtoMuNu_minus/WtoMuNu_minus.root PomwigWtoMuNu_minus/res/*.root
hadd WtoENu_plus/WtoENu_plus.root PomwigWtoENu_plus/res/*.root
hadd WtoMuNu_plus/WtoMuNu_plus.root PomwigWtoMuNu_plus/res/*.root
hadd ZtoMuMu_plus/ZtoMuMu_plus.root PomwigZtoMuMu_plus/res/*.root
hadd ZtoMuMu_plus/ZtoMuMu_plus.root PomwigZtoMuMu_plus/res/*.root
hadd ZtoMuMu_minus/ZtoMuMu_minus.root PomwigZtoMuMu_minus/res/*.root
hadd ZtoEE_minus/ZtoEE_minus.root PomwigZtoEE_minus/res/*.root
hadd ZtoEE_plus/ZtoEE_plus.root PomwigZtoEE_plus/res/*.root
