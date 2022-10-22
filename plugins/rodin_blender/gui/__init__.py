import bpy

import rodin_blender.gui.topbar as topbar

def register():
    # Top-bar menu
    bpy.utils.register_class(topbar.TOPBAR_MT_Rodin)
    bpy.utils.register_class(topbar.TOPBAR_MT_RodinExport)
    bpy.types.TOPBAR_MT_editor_menus.append(topbar.TOPBAR_MT_Rodin.menu_draw)

def unregister():
    # Top-bar menu
    bpy.types.TOPBAR_MT_editor_menus.remove(topbar.TOPBAR_MT_Rodin.menu_draw)
    bpy.utils.unregister_class(topbar.TOPBAR_MT_RodinExport)
    bpy.utils.unregister_class(topbar.TOPBAR_MT_Rodin)
