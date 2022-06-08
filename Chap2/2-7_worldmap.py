### Saizensen, Fig 2.7: drawing world map
### June 2022 / Masakiyo Kitazawa

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

plt.figure(linewidth=5)
ax = plt.axes(projection=ccrs.PlateCarree() )
#ax = plt.axes(projection=ccrs.Mollweide(),lw=0)
ax.set_extent([-170,170,-50,80], ccrs.PlateCarree() )

land_50m  = cfeature.NaturalEarthFeature('physical', 'land', '50m', 
                                         edgecolor='none',
                                         facecolor=cfeature.COLORS['land'])
ax.add_feature(land_50m)

countries_50m  = cfeature.NaturalEarthFeature('cultural', 'admin_0_countries', '50m', 
                                              edgecolor='gray',
                                              facecolor='none')
ax.add_feature( countries_50m , lw=0.1 )
ax.coastlines( resolution='50m' , lw=0.3 )
plt.savefig( "worldmap.pdf" )
