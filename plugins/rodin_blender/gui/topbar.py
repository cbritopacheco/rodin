import bpy

class TOPBAR_MT_Rodin(bpy.types.Menu):
    bl_label = "Rodin"

    def draw(self, context):
        layout = self.layout
        layout.menu('TOPBAR_MT_RodinExport')

    def menu_draw(self, context):
        self.layout.menu('TOPBAR_MT_Rodin')

class TOPBAR_MT_RodinExport(bpy.types.Menu):
    bl_label = 'Export'

    def draw(self, context):
        layout = self.layout
        layout.operator("mesh.primitive_cube_add")
